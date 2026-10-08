/**
 * @file temperature.cpp
 * @brief Implements region expansion and scalar temperature constraints.
 *
 * A prescribed temperature is a scalar thermal condition. The semantic target
 * may already be nodal or may require expansion through compiled element or
 * surface connectivity. Shared nodes are deduplicated before one equation
 *
 *     T_i = T_bar
 *
 * is appended per physical node in ascending node-ID order, independently
 * of the region traversal order.
 *
 * The implementation does not assemble heat flux, conductivity or convection
 * terms. It only constructs the algebraic primary-variable prescriptions that a
 * thermal analysis applies to its scalar temperature system.
 *
 * @see Temperature
 * @see Condition
 * @see constraint::Equation
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "temperature.h"

#include "../../core/config.h"
#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/geometry/surface/surface_interface.h"
#include "../../model/model_data.h"

#include <atomic>
#include <cmath>
#include <exception>
#include <sstream>
#include <vector>

namespace fem::bc {

/**
 * Resolves the selected region and prescribes one scalar thermal DOF per node.
 *
 * A node region already identifies the constrained degrees of freedom, whereas
 * surface and element regions must first be expanded through their compiled
 * nodal connectivity. Shared nodes are marked only once, even if several
 * selected surfaces or elements refer to them.
 *
 * Each unique target node i contributes the equation
 *
 *     1 * T_i = temperature_.
 *
 * The scalar thermal unknown occupies nodal component zero. Equations are
 * appended in ascending node-ID order, independent of parallel traversal
 * order. The stored temperature is absolute and is not interpolated or
 * amplitude-scaled by this condition.
 *
 * No RHS load, stiffness contribution or system numbering is modified.
 *
 * @param model_data Compiled nodal, element and surface connectivity.
 * @param rhs Scalar thermal RHS left unchanged.
 * @param equations Receives one prescribed-temperature row per unique node.
 * @param system_dof_ids Thermal DOF numbering left unchanged.
 * @param matrix Sparse matrix triplets left unchanged.
 * @param time Unused by this absolute-temperature prescription.
 * @param ignore_amplitude Unused by this condition.
 * @param step_progress Unused; no temperature history is interpolated.
 */
void Temperature::apply(
    model::ModelData&      model_data,
    model::Field&,
    constraint::Equations& equations,
    const SystemDofIds&,
    TripletList&,
    Precision,
    bool,
    Precision
) {
    // Exactly one semantic target is permitted. The resulting node indices
    // must refer to the compiled geometry used by the thermal system.
    const int active_regions = static_cast<int>(node_region_    != nullptr)
                             + static_cast<int>(surface_region_ != nullptr)
                             + static_cast<int>(element_region_ != nullptr);

    logging::error(active_regions == 1,
        "TEMPERATURE: exactly one node, surface or element region must be configured");
    logging::error(std::isfinite(temperature_),
        "TEMPERATURE: prescribed temperature must be finite");
    logging::error(model_data.positions != nullptr,
        "TEMPERATURE: model positions are not initialized");

    const Index node_count = model_data.positions->rows;

    // Atomic flags allow parallel region expansion even when multiple selected
    // elements or surfaces share the same physical node. Initialize all flags
    // before the parallel traversal so marking remains idempotent.
    std::vector<std::atomic_bool> selected(static_cast<std::size_t>(node_count));
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 4096) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index node = 0; node < node_count; ++node) {
        selected[static_cast<std::size_t>(node)].store(
            false,
            std::memory_order_relaxed
        );
    }

    std::atomic_bool   failed{false};
    std::exception_ptr failure = nullptr;

    auto mark_node = [&](ID node_id) {
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
            "TEMPERATURE: node ", node_id, " is outside the compiled node domain");

        selected[static_cast<std::size_t>(node_id)].store(
            true,
            std::memory_order_relaxed
        );
    };

    // Direct node selections need no topology expansion.
    if (node_region_) {
        const auto& node_ids = node_region_->data();
        const Index count = static_cast<Index>(node_ids.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static, 4096) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
        for (Index index = 0; index < count; ++index) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }

            try {
                mark_node(node_ids[static_cast<std::size_t>(index)]);
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
    }

    // Resolve element selections to their global nodal connectivity.
    if (element_region_) {
        const auto& element_ids = element_region_->data();
        const Index count = static_cast<Index>(element_ids.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
        for (Index index = 0; index < count; ++index) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }

            try {
                const ID element_id = element_ids[static_cast<std::size_t>(index)];
                logging::error(element_id >= 0
                            && static_cast<Index>(element_id) < static_cast<Index>(model_data.elements.size()),
                    "TEMPERATURE: element ", element_id,
                    " is outside the compiled element domain");

                const auto& element = model_data.elements[static_cast<std::size_t>(element_id)];
                logging::error(element != nullptr,
                    "TEMPERATURE: element ", element_id, " is not initialized");

                for (ID node_id : *element) {
                    mark_node(node_id);
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
    }

    // Resolve surface selections to their global nodal connectivity.
    if (surface_region_) {
        const auto& surface_ids = surface_region_->data();
        const Index count = static_cast<Index>(surface_ids.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
        for (Index index = 0; index < count; ++index) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }

            try {
                const ID surface_id = surface_ids[static_cast<std::size_t>(index)];
                logging::error(surface_id >= 0
                            && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
                    "TEMPERATURE: surface ", surface_id,
                    " is outside the compiled surface domain");

                const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
                logging::error(surface != nullptr,
                    "TEMPERATURE: surface ", surface_id, " is not initialized");

                for (ID node_id : *surface) {
                    mark_node(node_id);
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
    }

    if (failure) {
        std::rethrow_exception(failure);
    }

    // Count each selected node once before reserving the equation list.
    std::size_t selected_count = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:selected_count) schedule(static, 4096) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index node = 0; node < node_count; ++node) {
        if (selected[static_cast<std::size_t>(node)].load(std::memory_order_relaxed)) {
            ++selected_count;
        }
    }

    equations.reserve(equations.size() + selected_count);

    // Emit the prescribed temperature equation [1] * T_i = T_bar.
    // Serial emission guarantees stable node-ID ordering despite parallel
    // selection, and avoids synchronization on the equation collection.
    for (Index node = 0; node < node_count; ++node) {
        if (!selected[static_cast<std::size_t>(node)].load(std::memory_order_relaxed)) {
            continue;
        }

        const constraint::EquationEntry entry{
            static_cast<ID>(node),
            Dim(0),
            Precision(1)
        };
        equations.emplace_back(
            std::initializer_list<constraint::EquationEntry>{entry},
            temperature_
        );
    }
}


/**
 * Builds a compact diagnostic representation of the prescribed temperature.
 *
 * The output reports the semantic region without expanding its connectivity and
 * includes the absolute temperature stored by the condition.
 *
 * @return Human-readable temperature-condition description.
 */
std::string Temperature::str() const {
    // Report the semantic target without traversing or expanding connectivity.
    std::string target = "(unknown)";

    if (node_region_) {
        target = "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
    } else if (surface_region_) {
        target = "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
    } else if (element_region_) {
        target = "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
    }

    std::ostringstream os;
    os << "TEMPERATURE: target=" << target << ", value=" << temperature_;
    return os.str();
}

} // namespace fem::bc
