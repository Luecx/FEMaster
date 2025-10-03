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
#include "../core/logging.h"

#include <iomanip>
#include <sstream>

namespace fem::constraint {

namespace {
constexpr double kEps = 1e-8; // hardcoded epsilon as requested
}

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
    (void) model_data;
    Equations equations{};

    constexpr double kEps = 1e-8; // hardcoded epsilon threshold

    const char axis_names[3] = {'X', 'Y', 'Z'};

    for (int i = 0; i < 6; ++i) {
        if (!constrained_dofs_(i)) continue;

        // local basis for this DOF (Tx,Ty,Tz,Rx,Ry,Rz)
        const Vec3 local_dof_direction = (i < 3) ? Vec3::Unit(i) : Vec3::Unit(i - 3);

        // project to global
        const Vec3 g = coordinate_system_->to_global(local_dof_direction);

        std::vector<EquationEntry> entries;
        const int block = (i < 3) ? 0 : 1; // 0: translations, 1: rotations

        for (int j = 0; j < 3; ++j) {
            const double c = g(j);
            if (std::abs(c) <= kEps) continue; // drop tiny components

            const Dim dof = 3 * block + j;

            const auto dof_id_n1 = system_nodal_dofs(node_1_id_, dof);
            const auto dof_id_n2 = system_nodal_dofs(node_2_id_, dof);
            if (dof_id_n1 < 0 || dof_id_n2 < 0) {
                // silently skip this component; it's not active in the system
                continue;
            }

            entries.emplace_back(EquationEntry{node_1_id_, dof, (Precision) c});
            entries.emplace_back(EquationEntry{node_2_id_, dof, (Precision)-c});
        }

        // Only push if something meaningful remains
        if (!entries.empty()) {
            equations.emplace_back(Equation{std::move(entries)});
        }

        logging::warning(entries.empty(),
                         "Connector (", node_1_id_, " <-> ", node_2_id_,
                         ") generated no equation for local ",
                         (i < 3 ? "T" : "R"), axis_names[(i < 3) ? i : i - 3],
                         " (", g(0), ", ", g(1), ", ", g(2), "). ",
                         "in coordinate system '", coordinate_system_->name, "'.");
    }

    return equations;
}


/******************************************************************************
 * @brief Returns the coupled DOFs for the two nodes based on the connector.
 *
 * @return DOFs that are constrained by the connector.
 ******************************************************************************/
Dofs Connector::dofs() const {
    // Output: which GLOBAL DOFs (Tx,Ty,Tz,Rx,Ry,Rz) the connector actually impacts
    Dofs impacted;
    impacted.setZero();

    // For each local DOF we intend to constrain, transform its local basis direction to global.
    for (int i = 0; i < 6; ++i) {
        if (!constrained_dofs_(i)) continue; // skip unconstrained local directions

        // Build local unit direction for this DOF (Tx,Ty,Tz,Rx,Ry,Rz)
        // For rotations we use the same axis unit vectors; the coordinate system should map them.
        Vec3 local_dir = (i < 3) ? Vec3::Unit(i) : Vec3::Unit(i - 3);

        // Map local axis to global
        Vec3 g = coordinate_system_->to_global(local_dir);

        // Decide which GLOBAL components are significantly present
        const int block = (i < 3) ? 0 : 1; // 0 = translations, 1 = rotations
        for (int j = 0; j < 3; ++j)
        {
            if (std::abs(g(j)) > kEps)
            {
                const int global_dof = 3 * block + j; // (Tx,Ty,Tz,Rx,Ry,Rz)
                impacted(global_dof) = true;
            }
        }
    }

    return impacted;
}


}    // namespace fem::constraint
