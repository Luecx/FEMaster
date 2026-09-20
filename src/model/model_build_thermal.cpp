/**
 * @file model_build_thermal.cpp
 * @brief Implements global thermal DOF, operator, load and constraint assembly.
 *
 * The thermal system uses one scalar primary variable per active node. Volume
 * conduction is assembled from the conductivity matrices supplied by
 * `ThermalElement`, while the selected `ThermalCollector` objects contribute
 * prescribed heat flow, mixed boundary operators and temperature constraints.
 *
 * For a stationary conduction problem the resulting unconstrained system is
 *
 *     (K_T + K_b) T = q_b,
 *
 * where
 *
 *     K_T = sum_e integral_Omega_e grad(N)^T k grad(N) dOmega
 *
 * is the material conductivity operator, `K_b` contains unknown-dependent mixed
 * boundary terms such as convection, and `q_b` contains prescribed thermal
 * boundary sources. Essential temperatures remain separate equations `C T = d`
 * and are applied later by the constraint transformer.
 *
 * @see ThermalElement
 * @see bc::ThermalCollector
 * @see Model::build_thermal_dof_index_matrix
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#include "../core/config.h"
#include "../core/logging.h"
#include "../mattools/assemble.h"
#include "../mattools/numerate_dofs.h"
#include "element/element_thermal.h"
#include "model.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace fem::model {

/**
 * Enumerates the scalar temperature DOFs required by thermal elements.
 *
 * Every node referenced by at least one `ThermalElement` receives exactly one
 * active temperature component. Nodes referenced only by non-thermal elements
 * remain inactive. The resulting boolean mask is converted to contiguous
 * zero-based global thermal equation identifiers.
 *
 * @return Node-by-one matrix of global thermal system DOF identifiers.
 */
SystemDofIds Model::build_thermal_dof_index_matrix() {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    const Index node_count    = _data->positions->rows;
    const Index element_count = static_cast<Index>(_data->elements.size());

    // Several neighboring thermal elements may activate the same node
    // concurrently. Atomic flags make this idempotent write safe without one
    // complete nodal mask per worker.
    std::vector<std::atomic_bool> active_nodes(static_cast<std::size_t>(node_count));
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 4096) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index node = 0; node < node_count; ++node) {
        active_nodes[static_cast<std::size_t>(node)].store(
            false,
            std::memory_order_relaxed
        );
    }

    std::atomic_bool   failed{false};
    std::exception_ptr failure = nullptr;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        if (failed.load(std::memory_order_relaxed)) {
            continue;
        }

        try {
            const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
            if (element == nullptr || element->as<ThermalElement>() == nullptr) {
                continue;
            }

            for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
                const ID node_id = element->nodes()[local_node];
                logging::error(node_id >= 0 && static_cast<Index>(node_id) < node_count,
                    "Model: thermal element ", element->elem_id,
                    " references node ", node_id, " outside the compiled node domain");

                active_nodes[static_cast<std::size_t>(node_id)].store(
                    true,
                    std::memory_order_relaxed
                );
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

    SystemDofs mask{node_count, 1};
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 4096) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index node = 0; node < node_count; ++node) {
        mask(node, 0) = active_nodes[static_cast<std::size_t>(node)].load(
            std::memory_order_relaxed
        );
    }

    // Contiguous global equation numbering remains a deterministic prefix-style
    // operation after the parallel activation phase.
    return mattools::numerate_dofs(mask);
}

/**
 * Assembles the global material conductivity matrix.
 *
 * Every thermally participating element supplies
 *
 *     K_T^e = integral_Omega_e grad(N)^T k grad(N) dOmega.
 *
 * The common sparse matrix assembler maps the element-local scalar temperature
 * rows and columns through `system_dof_ids` and sums all overlapping element
 * contributions. Non-thermal elements are excluded before assembly so the local
 * callback always returns a valid square thermal element matrix.
 *
 * @param system_dof_ids Scalar node-to-system thermal equation mapping.
 * @return Assembled global conductivity matrix `K_T`.
 */
SparseMatrix Model::build_thermal_conductivity_matrix(
    const SystemDofIds& system_dof_ids
) {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");
    logging::error(system_dof_ids.rows() == _data->positions->rows,
        "Model: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "Model: thermal DOF map must contain exactly one component");

    // Select thermally participating elements in parallel. The indexed temporary
    // vector avoids synchronized push_back operations; compaction is linear and
    // cheap compared with element integration and sparse assembly.
    const Index element_count = static_cast<Index>(_data->elements.size());
    std::vector<ElementPtr> thermal_elements(static_cast<std::size_t>(element_count));

#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 2048) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        const auto& element = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (element != nullptr && element->as<ThermalElement>() != nullptr) {
            thermal_elements[static_cast<std::size_t>(elem_idx)] = element;
        }
    }

    thermal_elements.erase(
        std::remove(thermal_elements.begin(), thermal_elements.end(), nullptr),
        thermal_elements.end()
    );

    // The common assembler already parallelizes local conductivity evaluation,
    // triplet generation and sparse k-way merging when OpenMP is enabled.
    return mattools::assemble_matrix(
        thermal_elements,
        system_dof_ids,
        [](const ElementPtr& element, Precision* buffer) {
            auto* thermal = element->as<ThermalElement>();
            logging::error(thermal != nullptr,
                "Model: non-thermal element reached thermal conductivity assembly");
            return thermal->conductivity(buffer);
        }
    );
}

/**
 * Assembles the prescribed thermal boundary contribution into a scalar nodal RHS.
 *
 * Each selected `ThermalCollector` superimposes its load-like conditions into
 * one one-component NODE field. Pure Neumann heat flux contributes directly, and
 * Mixed conditions contribute only their prescribed source part during this
 * pass. Unknown-dependent Mixed terms are assembled separately by
 * `build_thermal_boundary_matrix()`.
 *
 * The resulting field is intentionally kept in nodal form so the later thermal
 * loadcase can reduce it through the active scalar DOF map using the same common
 * matrix-to-vector utility as structural analyses.
 *
 * @param thermal_sets Names of thermal collectors participating in the analysis.
 * @param time Analysis time forwarded to amplitude-dependent boundary conditions.
 * @return Scalar nodal thermal RHS field.
 */
Field Model::build_thermal_load_matrix(
    const std::vector<std::string>& thermal_sets,
    Precision time
) {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    Field rhs{
        "THERMAL_LOAD",
        FieldDomain::NODE,
        _data->field_rows(FieldDomain::NODE),
        1
    };
    rhs.set_zero();

    // Superimpose all selected thermal boundary source terms in user-provided
    // collector order.
    for (const std::string& name : thermal_sets) {
        logging::error(_data->thermal_cols.has(name),
            "Model: thermal collector ", name, " does not exist");

        const auto collector = _data->thermal_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal collector ", name, " is not initialized");

        collector->apply_rhs(*_data, rhs, time);
    }

    return rhs;
}

/**
 * Assembles the unknown-dependent operator contribution of mixed thermal BCs.
 *
 * Mixed conditions such as convection contribute a boundary matrix of the form
 *
 *     K_b^e = integral_Gamma_e h N N^T dGamma.
 *
 * Each selected collector appends sparse triplets in the active scalar thermal
 * numbering. Duplicate entries generated by neighboring boundary faces are
 * summed by Eigen while the final sparse matrix is constructed.
 *
 * @param system_dof_ids Scalar node-to-system thermal equation mapping.
 * @param thermal_sets Names of thermal collectors participating in the analysis.
 * @param time Analysis time forwarded to amplitude-dependent boundary conditions.
 * @return Assembled mixed thermal boundary operator `K_b`.
 */
SparseMatrix Model::build_thermal_boundary_matrix(
    const SystemDofIds&             system_dof_ids,
    const std::vector<std::string>& thermal_sets,
    Precision                       time
) {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");
    logging::error(system_dof_ids.rows() == _data->positions->rows,
        "Model: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "Model: thermal DOF map must contain exactly one component");

    const int system_size = system_dof_ids.size() == 0
        ? 0
        : system_dof_ids.maxCoeff() + 1;

    TripletList triplets;

    // Mixed conditions append their local boundary operators directly in active
    // global thermal equation numbering.
    for (const std::string& name : thermal_sets) {
        logging::error(_data->thermal_cols.has(name),
            "Model: thermal collector ", name, " does not exist");

        const auto collector = _data->thermal_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal collector ", name, " is not initialized");

        collector->apply_matrix(*_data, system_dof_ids, triplets, time);
    }

    SparseMatrix matrix(system_size, system_size);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    return matrix;
}

/**
 * Collects prescribed-temperature equations from selected thermal collectors.
 *
 * Thermal essential conditions are kept distinct from structural supports and
 * kinematic constraints. For the current thermal boundary-condition set each
 * `Temperature` object contributes rows
 *
 *     T_i = T_bar.
 *
 * The returned equations are not transformed or reduced here. A thermal
 * loadcase can pass them directly to `ConstraintTransformer` after the complete
 * material and mixed boundary operators have been assembled.
 *
 * @param thermal_sets Names of thermal collectors participating in the analysis.
 * @return Concatenated scalar thermal Dirichlet equations.
 */
constraint::Equations Model::collect_thermal_constraints(
    const std::vector<std::string>& thermal_sets
) {
    constraint::Equations equations;

    for (const std::string& name : thermal_sets) {
        logging::error(_data->thermal_cols.has(name),
            "Model: thermal collector ", name, " does not exist");

        const auto collector = _data->thermal_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal collector ", name, " is not initialized");

        auto collector_equations = collector->get_equations(*_data);
        equations.reserve(equations.size() + collector_equations.size());

        for (auto& equation : collector_equations) {
            equations.push_back(std::move(equation));
        }
    }

    return equations;
}

} // namespace fem::model
