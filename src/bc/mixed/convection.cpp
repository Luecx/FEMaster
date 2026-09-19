/**
 * @file convection.cpp
 * @brief Implements linear thermal convection RHS and boundary-operator assembly.
 *
 * Newton cooling separates naturally into a prescribed ambient source
 *
 *     q_h = integral_Gamma h T_inf N^T dGamma
 *
 * and an unknown-dependent symmetric boundary operator
 *
 *     K_h = integral_Gamma h N N^T dGamma.
 *
 * Both terms are integrated in the reference configuration. Shape-product
 * integration uses a dedicated higher-order surface quadrature because
 * `N_i N_j` has higher polynomial order than ordinary surface loading.
 *
 * @see Convection
 * @see model::SurfaceInterface::integrate_scalar_shape_matrix
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#include "convection.h"

#include "../../core/config.h"
#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <sstream>
#include <vector>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem::bc {

/**
 * Evaluates the effective film coefficient used by both convection contributions.
 *
 * The nominal film coefficient must be finite and non-negative. If an amplitude
 * is active, it scales `h`; the resulting effective coefficient is validated
 * again because an arbitrary amplitude could otherwise reverse the dissipative
 * boundary operator or introduce a non-finite value.
 *
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Use the nominal coefficient directly when true.
 * @return Validated effective film coefficient.
 */
Precision Convection::effective_film_coefficient(Precision time, bool ignore_amplitude) const {
    // Validate the physical parameters before applying temporal scaling
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");
    logging::error(std::isfinite(ambient_temperature_),
        "CONVECTION: ambient temperature must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision h = film_coefficient_ * scale;

    logging::error(std::isfinite(h) && h >= Precision(0),
        "CONVECTION: effective film coefficient must be finite and non-negative");

    return h;
}

/**
 * Assembles the prescribed ambient source into the scalar thermal RHS.
 *
 * The source term of Newton cooling is
 *
 *     q_h = integral_Gamma h T_inf N^T dGamma.
 *
 * It is independent of the unknown temperature field and can therefore be
 * assembled through the ordinary load-like RHS interface.
 *
 * @param model_data Compiled surface topology and reference geometry.
 * @param rhs Scalar nodal thermal RHS receiving the ambient source.
 * @param time Analysis time used for optional amplitude evaluation.
 * @param ignore_amplitude Apply the nominal film coefficient when true.
 */
void Convection::apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude) {
    // Validate the target and scalar thermal assembly context before entering the
    // parallel surface traversal.
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "CONVECTION: target field must be a NODE field with exactly one component");

    // A vanishing film coefficient removes both convection contributions exactly.
    const Precision h = effective_film_coefficient(time, ignore_amplitude);
    if (h == Precision(0)) {
        return;
    }

    const Precision source = h * ambient_temperature_;
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
    // Integrate each surface independently and scatter only its local nodal
    // vector. Atomic scalar additions preserve shared-edge contributions without
    // one full nodal RHS allocation per worker.
    for (Index surface_index = 0; surface_index < surface_count; ++surface_index) {
        if (failed.load(std::memory_order_relaxed)) {
            continue;
        }

        try {
            const ID surface_id = surface_ids[static_cast<std::size_t>(surface_index)];
            logging::error(surface_id >= 0
                        && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
                "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

            const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
            logging::error(surface != nullptr,
                "CONVECTION: surface ", surface_id, " is not initialized");

            const DynamicVector local = surface->integrate_scalar_shape_vector(
                *model_data.positions_reference,
                [source](const Vec3&) -> Precision { return source; }
            );

            logging::error(local.size() == surface->n_nodes && local.allFinite(),
                "CONVECTION: local ambient source is invalid on surface ", surface_id);

            for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
                const ID node_id = surface->nodes()[local_node];
                logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
                    "CONVECTION: surface references node ", node_id,
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
 * Assembles the temperature-dependent convection boundary operator.
 *
 * For each selected surface the local matrix
 *
 *     K_h^e = integral_Gamma_e h N N^T dGamma
 *
 * is integrated in connectivity ordering. The scalar thermal DOF map converts
 * every surface node to an active system row and column before non-zero entries
 * are appended to the global sparse triplet list.
 *
 * @param model_data Compiled surface topology and reference geometry.
 * @param system_dof_ids Scalar node-to-active-temperature equation mapping.
 * @param matrix Sparse triplet list receiving convection operator entries.
 * @param time Analysis time used for optional amplitude evaluation.
 * @param ignore_amplitude Apply the nominal film coefficient when true.
 */
void Convection::apply_matrix(model::ModelData&   model_data,
                              const SystemDofIds& system_dof_ids,
                              TripletList&        matrix,
                              Precision           time,
                              bool                ignore_amplitude) {
    // Validate the surface target and scalar thermal system numbering before
    // parallel local-operator evaluation.
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(system_dof_ids.rows() == model_data.positions_reference->rows,
        "CONVECTION: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "CONVECTION: thermal DOF map must contain exactly one component");

    const Precision h = effective_film_coefficient(time, ignore_amplitude);
    if (h == Precision(0)) {
        return;
    }

    const auto& positions     = *model_data.positions_reference;
    const auto& surface_ids   = region_->data();
    const Index surface_count = static_cast<Index>(surface_ids.size());

#ifdef _OPENMP
    const int worker_count = std::max(
        1,
        std::min(global_config.max_threads, static_cast<int>(std::max<Index>(surface_count, 1)))
    );
#else
    const int worker_count = 1;
#endif

    // Each OpenMP worker owns its triplet buffer, eliminating synchronization
    // from the quadratic local matrix scatter. Duplicate global entries are
    // intentionally preserved for the final SparseMatrix construction.
    std::vector<TripletList> thread_triplets(static_cast<std::size_t>(worker_count));

    std::atomic_bool   failed{false};
    std::exception_ptr failure = nullptr;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 128) num_threads(worker_count) if(worker_count > 1)
#endif
    for (Index surface_index = 0; surface_index < surface_count; ++surface_index) {
        if (failed.load(std::memory_order_relaxed)) {
            continue;
        }

        try {
#ifdef _OPENMP
            const int worker = omp_get_thread_num();
#else
            const int worker = 0;
#endif
            auto& local_triplets = thread_triplets[static_cast<std::size_t>(worker)];

            const ID surface_id = surface_ids[static_cast<std::size_t>(surface_index)];
            logging::error(surface_id >= 0
                        && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
                "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

            const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
            logging::error(surface != nullptr,
                "CONVECTION: surface ", surface_id, " is not initialized");

            const DynamicMatrix local = surface->integrate_scalar_shape_matrix(
                positions,
                [h](const Vec3&) -> Precision { return h; }
            );

            logging::error(local.rows() == surface->n_nodes && local.cols() == surface->n_nodes,
                "CONVECTION: local boundary matrix does not match surface connectivity");
            logging::error(local.allFinite(),
                "CONVECTION: local boundary matrix contains NaN or Inf");

            for (Index i = 0; i < surface->n_nodes; ++i) {
                const ID node_i = surface->nodes()[i];
                logging::error(node_i >= 0
                            && static_cast<Eigen::Index>(node_i) < system_dof_ids.rows(),
                    "CONVECTION: surface references invalid node ", node_i);

                const int row = system_dof_ids(static_cast<Eigen::Index>(node_i), 0);
                logging::error(row >= 0,
                    "CONVECTION: surface references thermally inactive node ", node_i);

                for (Index j = 0; j < surface->n_nodes; ++j) {
                    const ID node_j = surface->nodes()[j];
                    logging::error(node_j >= 0
                                && static_cast<Eigen::Index>(node_j) < system_dof_ids.rows(),
                        "CONVECTION: surface references invalid node ", node_j);

                    const int col = system_dof_ids(static_cast<Eigen::Index>(node_j), 0);
                    logging::error(col >= 0,
                        "CONVECTION: surface references thermally inactive node ", node_j);

                    const Precision value = local(
                        static_cast<Eigen::Index>(i),
                        static_cast<Eigen::Index>(j)
                    );

                    if (value != Precision(0)) {
                        local_triplets.emplace_back(row, col, value);
                    }
                }
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

    // Concatenate the independent worker buffers once. Sparse construction later
    // combines duplicate row/column entries from adjacent surfaces.
    std::size_t additional = 0;
    for (const auto& local : thread_triplets) {
        additional += local.size();
    }
    matrix.reserve(matrix.size() + additional);

    for (const auto& local : thread_triplets) {
        matrix.insert(matrix.end(), local.begin(), local.end());
    }
}

/**
 * Builds the diagnostic representation of the convection condition.
 *
 * @return Human-readable target region, nominal film coefficient, ambient
 *         temperature and optional amplitude.
 */
std::string Convection::str() const {
    std::ostringstream os;

    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", ambient=" << ambient_temperature_;

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
