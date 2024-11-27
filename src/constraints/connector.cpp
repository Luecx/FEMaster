/******************************************************************************
 * @file connector.cpp
 * @brief Implements the Connector class methods for applying connector constraints between nodes in FEM.
 *
 * This file contains the implementation of the `get_equations` method, which generates the
 * constraint equations for linking two nodes with specific degrees of freedom (DOFs).
 *
 * @see connector.h
 * @date 08.10.2024
 ******************************************************************************/

#include "connector.h"

#include "../model/model_data.h"

namespace fem::constraint {

/******************************************************************************
 * @brief Constructor for the Connector class.
 *
 * Initializes a Connector object with the specified nodes, coordinate system, and connector type.
 * The constrained degrees of freedom (DOFs) are derived from the connector type using bitmasking.
 ******************************************************************************/
Connector::Connector(ID                              node_1_id,
                     ID                              node_2_id,
                     fem::cos::CoordinateSystem::Ptr coordinate_system,
                     ConnectorType                   connector_type)
    : coordinate_system_(coordinate_system)
    , connector_type_(connector_type)
    , constrained_dofs_({!!(connector_type_ & 0b100000),    // Translation-X
                         !!(connector_type_ & 0b010000),    // Translation-Y
                         !!(connector_type_ & 0b001000),    // Translation-Z
                         !!(connector_type_ & 0b000100),    // Rotation-X
                         !!(connector_type_ & 0b000010),    // Rotation-Y
                         !!(connector_type_ & 0b000001)})
    ,    // Rotation-Z
    node_1_id_(node_1_id)
    , node_2_id_(node_2_id) {}

/******************************************************************************
 * @brief Generates the constraint equations for the connector element.
 *
 * This method constructs the constraint equations for the two linked nodes based on the
 * constrained DOFs defined by the connector type. It transforms the local DOF directions
 * into the global system using the custom coordinate system.
 *
 * @param system_nodal_dofs The global system DOF IDs for each node.
 * @param node_coords The coordinates of all nodes (unused in this implementation).
 * @param row_offset The starting row for inserting the constraint equations.
 * @return TripletList A list of triplets representing the constraint equations.
 ******************************************************************************/
Equations Connector::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    (void) model_data;    // Unused parameter
    Equations equations {};

    // Loop through each degree of freedom (DOF)
    for (int i = 0; i < 6; i++) {
        if (!constrained_dofs_(i))
            continue;    // Skip unconstrained DOFs.

        // Get the global system DOF indices for the two nodes
        auto dof_1 = system_nodal_dofs(node_1_id_, i);
        auto dof_2 = system_nodal_dofs(node_2_id_, i);

        // Skip if any of the DOFs are invalid (-1 means unconstrained in the system)
        if (dof_1 < 0 || dof_2 < 0)
            continue;

        // Determine the local axis based on whether it's a translational or rotational DOF
        Vec3 local_dof_direction;
        if (i < 3) {
            // Translational DOFs (0: X, 1: Y, 2: Z)
            local_dof_direction = Vec3::Unit(i);
        } else {
            // Rotational DOFs (3: RX, 4: RY, 5: RZ)
            local_dof_direction = Vec3::Unit(i - 3);    // Adjust to correct local rotational axis
        }

        // Transform the local DOF direction to the global coordinate system
        Vec3                       global_dof_direction = coordinate_system_->to_global(local_dof_direction);

        std::vector<EquationEntry> entries;

        // Add the constraint equation to the triplet list
        for (int j = 0; j < 3; j++) {

            Dim  dof       = 3 * (i / 3) + j;

            auto dof_id_n1 = system_nodal_dofs(node_1_id_, dof);
            auto dof_id_n2 = system_nodal_dofs(node_2_id_, dof);

            if (dof_id_n1 < 0 || dof_id_n2 < 0) {
                logging::warning(false,
                                 "Invalid DOF indices for node ",
                                 node_1_id_,
                                 " and ",
                                 node_2_id_,
                                 ". Constraint may not work as intended.");
                continue;
            }

            EquationEntry en1 = EquationEntry {node_1_id_, dof, global_dof_direction(j)};
            EquationEntry en2 = EquationEntry {node_2_id_, dof, -global_dof_direction(j)};

            entries.push_back(en1);
            entries.push_back(en2);
        }
        equations.push_back(Equation {entries});
    }

    return equations;
}

}    // namespace fem::constraint
