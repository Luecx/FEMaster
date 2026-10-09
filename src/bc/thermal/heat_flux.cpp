/**
 * @file heat_flux.cpp
 * @brief Implements consistent prescribed thermal surface-flux assembly.
 *
 * The prescribed scalar heat flux is integrated in the reference configuration.
 * At every quadrature point the surface shape functions distribute the physical
 * heat-flow density to the connected nodes, producing
 *
 *     q_e = integral_Gamma_e N^T q_bar dGamma.
 *
 * Optional temporal amplitude scaling is evaluated once for the complete
 * condition before its target surfaces are traversed.
 *
 * @see HeatFlux
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#include "heat_flux.h"

#include "../../core/config.h"
#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the prescribed scalar surface heat flux into the thermal RHS.
 *
 * The flux is interpolated from its initial to target value using normalized
 * step progress, unless an assigned amplitude scales the target at physical
 * time. ignore_amplitude uses the nominal target without interpolation.
 *
 * For each selected surface Gamma_e, the consistent nodal heat-flow vector is
 *
 *     q_e = integral_Gamma_e N^T(x) * q dGamma
 *
 * and the quadrature implementation evaluates, for every surface node i,
 *
 *     q_i += sum_k N_i(x_k) * q * J_s(x_k) * w_k.
 *
 * The surface integration routine computes the shape functions N_i, physical
 * reference-area Jacobian J_s and quadrature weights w_k. This method scatters
 * the resulting local scalar contributions into the global nodal RHS. Adjacent
 * surfaces may share nodes, so OpenMP scattering uses atomic accumulation.
 *
 * @param model_data Compiled surfaces and reference nodal positions.
 * @param rhs Scalar nodal thermal RHS modified in place.
 * @param equations Constraint equations left unchanged.
 * @param system_dof_ids Global thermal DOF numbering left unchanged.
 * @param matrix Sparse matrix triplets left unchanged.
 * @param time Physical time used for amplitude evaluation.
 * @param ignore_amplitude Assemble the nominal target flux when true.
 * @param step_progress Normalized progress from initial to target flux.
 */
void HeatFlux::apply(
    model::ModelData&      model_data,
    model::Field&          rhs,
    constraint::Equations&,
    const SystemDofIds&,
    TripletList&,
    Precision              time,
    bool                   ignore_amplitude,
    Precision              step_progress
) {
    logging::error(region_ != nullptr,
        "HEATFLUX: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HEATFLUX: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "HEATFLUX: target field must be a NODE field with exactly one component");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    // Without an amplitude, q(lambda) = (1-lambda)*q_start + lambda*q_end.
    // An assigned amplitude instead evaluates q(time) = A(time)*q_end.
    Precision start_scale = Precision(1) - step_progress;
    Precision end_scale   = step_progress;

    if (ignore_amplitude) {
        start_scale = Precision(0);
        end_scale   = Precision(1);
    } else if (amplitude_) {
        start_scale = Precision(0);
        end_scale   = amplitude_->evaluate(time);
    }

    // An unspecified initial flux represents a zero starting contribution.
    const Precision start = std::isfinite(heat_flux_start_)
        ? heat_flux_start_ : Precision(0);
    const Precision flux = start_scale * start + end_scale * heat_flux_;

    logging::error(std::isfinite(flux),
        "HEATFLUX: effective heat flux must be finite");

    if (flux == Precision(0)) return;

    // Integrate each surface independently; accumulate shared nodal RHS entries
    // atomically when the surface traversal runs in parallel.
    const auto& surface_ids   = region_->data();
    const Index surface_count = static_cast<Index>(surface_ids.size());
    const Index node_count    = rhs.rows;
    Precision*  rhs_data      = rhs.data();

#ifdef _OPENMP
    const int worker_count = std::max(
        1,
        std::min(global_config.max_threads, static_cast<int>(std::max<Index>(surface_count, 1)))
    );
#else
    const int worker_count = 1;
#endif

    std::atomic_bool   failed{false};
    std::exception_ptr failure = nullptr;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 256) num_threads(worker_count) if(worker_count > 1)
#endif
    // Each worker integrates a local surface vector; atomic additions preserve
    // contributions from adjacent surfaces sharing the same thermal node.
    for (Index surface_index = 0; surface_index < surface_count; ++surface_index) {
        if (failed.load(std::memory_order_relaxed)) {
            continue;
        }

        try {
            const ID surface_id = surface_ids[static_cast<std::size_t>(surface_index)];
            logging::error(surface_id >= 0
                        && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
                "HEATFLUX: surface ", surface_id, " is outside the compiled surface domain");

            const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
            logging::error(surface != nullptr,
                "HEATFLUX: surface ", surface_id, " is not initialized");

            const DynamicVector local = surface->integrate_scalar_shape_vector(
                *model_data.positions_reference,
                [flux](const Vec3&) -> Precision { return flux; }
            );

            logging::error(static_cast<Index>(local.size()) == surface->n_nodes && local.allFinite(),
                "HEATFLUX: local surface load is invalid on surface ", surface_id);

            for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
                const ID node_id = surface->nodes()[local_node];
                logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
                    "HEATFLUX: surface references node ", node_id,
                    " outside the thermal RHS domain");

                const Precision contribution = local(local_node);
#ifdef _OPENMP
                #pragma omp atomic update
#endif
                rhs_data[static_cast<std::size_t>(node_id)] += contribution;
            }
        } catch (...) {
            failed.store(true, std::memory_order_relaxed);
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                if (!failure) {
                    failure = std::current_exception();
                }
            }
        }
    }

    if (failure) {
        std::rethrow_exception(failure);
    }
}


/**
 * Builds the diagnostic representation of the prescribed surface heat flux.
 *
 * @return Human-readable target region, nominal value and optional amplitude.
 */
std::string HeatFlux::str() const {
    std::ostringstream os;

    // Report the stored flux before temporal scaling or surface integration.
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), value=" << heat_flux_;

    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
