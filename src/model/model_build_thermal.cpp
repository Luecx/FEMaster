/**
 * @file model_build_thermal.cpp
 * @brief Implements global thermal DOF, operator, load and constraint assembly.
 *
 * The thermal system uses one scalar primary variable per active node. Volume
 * conduction is assembled from the conductivity matrices supplied by
 * `ThermalElement`. Direct thermal histories in ModelData::conditions and optional
 * selected ThermalCollector entries contribute prescribed heat flow, convection
 * operators and temperature constraints through Condition::apply.
 *
 * For a stationary conduction problem the resulting unconstrained system is
 *
 *     (K_T + K_b) T = q_b,
 *
 * where
 *
 *     K_T = sum_e integral_Omega_e grad(N)^T k grad(N) dOmega
 *
 * is the material conductivity operator, `K_b` contains unknown-dependent
 * convection terms, and `q_b` contains prescribed thermal boundary sources.
 * Prescribed temperatures remain separate equations `C T = d`
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
    logging::error(static_cast<Index>(system_dof_ids.rows()) == _data->positions->rows,
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
 * Assembles the current scalar thermal source field.
 *
 * Direct HEAT_FLUX and CONVECTION history and collector-activated thermal conditions
 * superimpose prescribed heat flow in a one-component NODE field. Convection
 * supplies q_h = integral_Gamma h T_inf N^T dGamma through the common apply
 * interface. Its simultaneous matrix contribution is discarded here because the
 * thermal solver obtains that operator from build_thermal_boundary_matrix().
 *
 * A valid thermal numbering is supplied when evaluating mixed conditions. This
 * keeps Condition::apply complete without RTTI, capability bases or special
 * collector dispatch. The returned nodal source is reduced to active equations
 * by the thermal analysis using the common matrix-to-vector utilities.
 *
 * @param time Analysis time supplied to amplitude-dependent conditions.
 * @return Scalar nodal thermal RHS for the current model state.
 */
Field Model::build_thermal_load_matrix(Precision time, Precision step_progress) {
    // Validate the scalar nodal domain before preparing common assembly outputs
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    Field rhs{"THERMAL_LOAD", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 1};
    rhs.set_zero();

    const SystemDofIds    system_dof_ids = build_thermal_dof_index_matrix();
    constraint::Equations equations{};
    TripletList           matrix{};

    // Heat flux and convection occupy separate active history families;
    // the parser has already selected any applicable collector definitions.
    for (const auto family : {bc::HEAT_FLUX, bc::CONVECTION}) {
        for (const auto& condition : _data->conditions.get(family)) {
            condition->apply(*_data, rhs, equations, system_dof_ids, matrix,
                             time, false, step_progress);
        }
    }

    return rhs;
}

/**
 * Assembles the current unknown-dependent thermal boundary operator.
 *
 * Active CONVECTION conditions append triplets for the boundary operator
 * K_h = integral_Gamma h N N^T dGamma in active scalar equation numbering.
 * Condition::apply also evaluates the ambient source into a temporary nodal RHS;
 * only the matrix is retained by this pass. The reader has already grouped
 * active thermal definitions into the appropriate ConditionManager families.
 *
 * Eigen sums overlapping entries while constructing the final sparse operator.
 * Temperature and heat-flux conditions are assembled in their own families.
 *
 * @param system_dof_ids Scalar node-to-system thermal equation mapping.
 * @param time Analysis time supplied to amplitude-dependent conditions.
 * @return Thermal boundary operator K_h for the current model state.
 */
SparseMatrix Model::build_thermal_boundary_matrix(const SystemDofIds& system_dof_ids, Precision time) {
    // Validate the scalar equation numbering and compiled nodal domain
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");
    logging::error(static_cast<Index>(system_dof_ids.rows()) == _data->positions->rows,
        "Model: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "Model: thermal DOF map must contain exactly one component");

    const int system_size = system_dof_ids.size() == 0 ? 0 : system_dof_ids.maxCoeff() + 1;
    Field rhs{"THERMAL_BOUNDARY_SOURCE", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 1};
    rhs.set_zero();

    constraint::Equations equations{};
    TripletList           triplets{};

    // Active convection contributes both source and film operator; retain K_h
    for (const auto& condition : _data->conditions.get(bc::CONVECTION)) {
        condition->apply(*_data, rhs, equations, system_dof_ids, triplets, time);
    }

    // Sum duplicate boundary entries into the global active thermal operator
    SparseMatrix matrix(system_size, system_size);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    return matrix;
}

/**
 * Collects current prescribed-temperature equations C T = d.
 *
 * Active TEMPERATURE conditions contribute scalar prescriptions T_i = T_bar.
 * These equations remain separate from
 * structural supports, MPCs and other structural kinematic constraints.
 *
 * Complete temporary RHS and numbering outputs preserve the common apply
 * interface. Contributions other than prescribed temperature equations are
 * discarded by this pass. The thermal analysis passes the collected equations
 * to its constraint transformer after assembling the material and film operators.
 *
 * @return Scalar temperature equations for the current model state.
 */
constraint::Equations Model::collect_thermal_constraints() {
    // Prepare all common outputs, retaining only prescribed-temperature rows
    Field rhs{"THERMAL_CONSTRAINT_SOURCE", FieldDomain::NODE, _data->field_rows(FieldDomain::NODE), 1};
    rhs.set_zero();

    const SystemDofIds    system_dof_ids = build_thermal_dof_index_matrix();
    constraint::Equations equations{};
    TripletList           matrix{};

    // Initial and analysis-local temperatures share the same persistent history
    for (const auto& condition : _data->conditions.get(bc::TEMPERATURE)) {
        condition->apply(*_data, rhs, equations, system_dof_ids, matrix, Precision(0), true);
    }

    return equations;
}

} // namespace fem::model
