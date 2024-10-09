//
// Created by f_eggers on 08.10.2024.
//

#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "../cos/coordinate_system.h"
#include "../core/types.h"

/**
 * @enum ConnectorType
 * @brief Enum for defining different types of connectors using bitwise integer representations.
 *
 * @details Each connector type is defined using a bitmask, where each bit represents a DOF:
 *          - Bit 0: Translation-X
 *          - Bit 1: Translation-Y
 *          - Bit 2: Translation-Z
 *          - Bit 3: Rotation-X
 *          - Bit 4: Rotation-Y
 *          - Bit 5: Rotation-Z
 */
enum ConnectorType : int {
    None  = 0b000000,     ///< No DOFs are constrained.
    Beam  = 0b111111,     ///< All 6 DOFs are constrained.
    Hinge = 0b111011      ///< All DOFs except rotation around X (bit 3 is 0).
};

/**
 * @class Connector
 * @brief Class for representing a connector element that links two nodes.
 *
 * @details This class constrains specific degrees of freedom (DOFs) between
 *          two nodes using a custom coordinate system. Each DOF can be independently
 *          constrained based on the user-defined `constrained_dofs_` vector.
 */
class Connector {
    fem::cos::CoordinateSystemPtr coordinate_system_; ///< Pointer to the custom coordinate system.
    ConnectorType connector_type_;          ///< Type of connector element.
    Dofs constrained_dofs_;                 ///< Vector indicating which DOFs are constrained.
    ID node_1_id_;                          ///< ID of the first node.
    ID node_2_id_;                          ///< ID of the second node.

public:
    /**
     * @brief Constructor for the Connector class.
     * @param node_1_id ID of the first node.
     * @param node_2_id ID of the second node.
     * @param coordinate_system Pointer to the custom coordinate system.
     * @param constrained_dofs Boolean vector indicating which DOFs are constrained.
     */
    Connector(ID node_1_id, ID node_2_id, fem::cos::CoordinateSystemPtr coordinate_system,  ConnectorType connector_type)
        : coordinate_system_(coordinate_system), connector_type_(connector_type),
        constrained_dofs_({
                !!(connector_type_ & 0b100000),
                !!(connector_type_ & 0b010000),
                !! (connector_type_ & 0b001000),
                !! (connector_type_ & 0b000100),
                !! (connector_type_ & 0b000010),
                !! (connector_type_ & 0b000001)}),
    node_1_id_(node_1_id), node_2_id_(node_2_id) {

    }

    virtual ~Connector() = default;

    /**
     * @brief Generate the constraint equations for the connector element.
     * @param system_nodal_dofs The global system DOF IDs for each node.
     * @param node_coords The coordinates of all nodes.
     * @param row_offset The starting row for the constraint equations.
     * @return A list of triplets representing the constraint equations.
     */
    TripletList get_equations(SystemDofIds& system_nodal_dofs, NodeData& /* node_coords */, int row_offset) {
        TripletList triplets{}; // Initialize the triplet list for constraint equations.
        int row = row_offset;   // Row index starts from the given offset.

        // Loop through each degree of freedom (DOF)
        for (int i = 0; i < 6; i++) {
            if (!constrained_dofs_(i))
                continue; // Skip unconstrained DOFs.

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
                local_dof_direction = Vec3::Unit(i - 3); // Adjust to correct local rotational axis
            }
            // could also use Vec3::Unit(i % 3);

            // Transform the local DOF direction to the global coordinate system
            Vec3 global_dof_direction = coordinate_system_->to_global(local_dof_direction);

            // Add the constraint equation to the triplet list
            for (int j = 0; j < 3; j++) {
                // i / 3 gives either the offset for the translational or rotational DOFs
                // and j is the index for the x, y, z components of the global DOF direction
                triplets.push_back(Triplet(row, i / 3 + j, global_dof_direction(j)));
                triplets.push_back(Triplet(row, i / 3 + j, -global_dof_direction(j)));
            }

            row++; // Move to the next row for the next constraint
        }

        return triplets; // Return the list of equations.
    }
};

#endif //CONNECTOR_H
