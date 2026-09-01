/**
 * @file model_thermal_load.cpp
 * @brief Implements model-side thermal boundary-condition assembly.
 *
 * The thermal load path constructs a one-component nodal heat-flow field and
 * dispatches selected load-collector entries according to their mathematical
 * boundary-condition type. Thermal Neumann conditions modify only the nodal RHS;
 * Robin conditions modify that same RHS and append sparse entries to the thermal
 * operator through the active scalar temperature DOF mapping.
 *
 * This model-level dispatch is the boundary between semantic load definitions
 * and analysis-specific algebraic assembly. Individual boundary conditions own
 * their local surface integration, while `Model` owns collector selection and
 * guarantees a common scalar thermal target layout.
 *
 * @see Model::build_thermal_load_matrix
 * @see bc::ThermalNeumann
 * @see bc::Robin
 * @see loadcase::SteadyStateThermal
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "model.h"

#include "../bc/neumann/neumann.h"
#include "../bc/robin/robin.h"
#include "../core/logging.h"

#include <memory>
#include <string>

namespace fem::model {

/**
 * Assembles selected thermal boundary conditions into RHS and operator terms.
 *
 * The function creates a scalar NODE field representing consistent nodal heat
 * flow. For every requested load collector, entries are classified dynamically:
 *
 * - `ThermalNeumann` entries assemble only into the scalar nodal RHS,
 * - `Robin` entries assemble into the same RHS and append operator triplets,
 * - structural or otherwise non-thermal entries are rejected explicitly.
 *
 * The supplied `system_dof_ids` must contain exactly one scalar component for
 * every compiled node. Robin conditions use that mapping to convert local
 * boundary matrices into global algebraic row and column indices. Neumann
 * conditions remain independent of the active algebraic ordering and write only
 * to the nodal field.
 *
 * @param load_sets Names of load collectors selected by the thermal loadcase.
 * @param system_dof_ids Scalar node-to-active-temperature system mapping.
 * @param lhs Sparse triplet list receiving Robin operator contributions.
 * @param time Analysis time passed to boundary-condition amplitudes.
 * @return One-component nodal field containing the assembled thermal RHS before
 *         reduction to active system ordering.
 */
Field Model::build_thermal_load_matrix(
    const std::vector<std::string>& load_sets,
    const SystemDofIds&             system_dof_ids,
    TripletList&                    lhs,
    Precision                       time
) {
    // Validate the compiled nodal domain and thermal system mapping before any
    // selected boundary condition receives assembly context.
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");
    logging::error(system_dof_ids.rows() == static_cast<Eigen::Index>(_data->positions->rows),
        "Model: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "Model: thermal DOF map must contain exactly one component per node");

    // Create the common scalar nodal heat-flow accumulator. Boundary conditions
    // scatter consistent contributions by global node identifier; reduction to
    // active algebraic ordering occurs later in the thermal loadcase.
    Field rhs{
        "THERMAL_LOAD",
        FieldDomain::NODE,
        _data->field_rows(FieldDomain::NODE),
        1
    };
    rhs.set_zero();

    // Apply the load collectors explicitly selected by the current thermal
    // loadcase in user-provided order.
    for (const std::string& name : load_sets) {
        // Resolve the named collector and fail clearly when the loadcase refers
        // to a missing or uninitialized set.
        logging::error(_data->load_cols.has(name),
            "Model: thermal load collector ", name, " does not exist");

        const auto collector = _data->load_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal load collector ", name, " is not initialized");

        // Dispatch heterogeneous load entries according to thermal boundary type.
        // Null entries are ignored consistently with the existing collector
        // behavior; every non-null incompatible type is reported as an error.
        for (const auto& load : collector->entries()) {
            if (!load) continue;

            if (auto neumann = std::dynamic_pointer_cast<bc::ThermalNeumann>(load)) {
                // Pure flux-type conditions need only the scalar nodal RHS field
                // and remain independent of active system indices.
                neumann->apply(*_data, rhs, time, false);
                continue;
            }

            if (auto robin = std::dynamic_pointer_cast<bc::Robin>(load)) {
                // Mixed conditions receive the same nodal RHS plus the active
                // scalar DOF map and global triplet list for their operator term.
                robin->apply(*_data, rhs, system_dof_ids, lhs, time, false);
                continue;
            }

            // A structural load inside a thermal load set would otherwise be
            // silently omitted and should therefore fail at the model boundary.
            logging::error(false,
                "Model: thermal load collector ", name,
                " contains non-thermal condition ", load->str());
        }
    }

    return rhs;
}

} // namespace fem::model
