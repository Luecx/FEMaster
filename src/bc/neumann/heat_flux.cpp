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
 * Integrates the prescribed heat flux over all selected surfaces.
 *
 * The target must be a one-component nodal field because thermal conduction has
 * one primary temperature DOF per active node. Surface geometry, quadrature,
 * shape-function weighting and nodal scattering are delegated to the generic
 * surface integration layer.
 *
 * @param model_data Compiled surface topology and reference geometry.
 * @param rhs Scalar nodal thermal RHS receiving consistent heat flow.
 * @param time Analysis time used for optional amplitude evaluation.
 * @param ignore_amplitude Apply the nominal unscaled value when true.
 */
void HeatFlux::apply(model::ModelData& model_data,
                     model::Field&     rhs,
                     Precision         time,
                     bool              ignore_amplitude) {
    // Validate the semantic target, reference geometry and scalar thermal field
    // before entering the parallel surface traversal.
    logging::error(region_ != nullptr,
        "HEATFLUX: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HEATFLUX: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "HEATFLUX: target field must be a NODE field with exactly one component");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    // Apply the common temporal amplitude once to the prescribed heat-flow density.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision value = heat_flux_ * scale;

    logging::error(std::isfinite(value),
        "HEATFLUX: effective heat flux must be finite");

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
    // Every worker integrates one complete surface into a local vector. Only the
    // final scalar nodal scatter is shared; atomic updates preserve contributions
    // from adjacent faces without allocating one full RHS field per thread.
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
                [value](const Vec3&) -> Precision { return value; }
            );

            logging::error(local.size() == surface->n_nodes && local.allFinite(),
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

    // Report the semantic target and unscaled physical value
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", value=" << heat_flux_;

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
