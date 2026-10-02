/**
 * @file mpc.h
 * @brief Declares simple multi-point constraints between two structural nodes.
 *
 * MPC objects represent the small set of Abaqus-style two-node kinematic
 * relations supported directly by FEMaster. They generate ordinary linear
 * constraint equations and therefore use the same constraint transformer as
 * supports, connectors, couplings and ties.
 *
 * BEAM uses rigid-link small-rotation kinematics with the second node acting as
 * the independent master. TIE equates all available generalized DOFs, while PIN
 * equates translations only.
 *
 * @see src/constraints/types/equation.h
 * @see src/constraints/types/coupling.h
 *
 * @author Finn Eggers
 * @date 01.10.2026
 */

#pragma once

#include "../../core/types_cls.h"
#include "../../core/types_eig.h"
#include "equation.h"

#include <cstdint>

namespace fem {
namespace model {
struct ModelData;
}
}

namespace fem {
namespace constraint {

/**
 * @brief Represents one two-node multi-point constraint.
 *
 * The first node is the dependent node and the second node is the independent
 * master, matching the Abaqus MPC ordering. Every supported type is converted
 * directly into homogeneous linear equations:
 *
 * - BEAM: rigid-link translations including the master rotation offset and
 *   equal available rotations,
 * - TIE: equal available translations and rotations,
 * - PIN: equal available translations.
 *
 * The class stores only the semantic MPC definition. Equation assembly,
 * reduction and multiplier handling remain responsibilities of the common
 * constraint subsystem.
 */
class Mpc {
public:
    enum class Type : std::uint8_t {
        Beam,
        Tie,
        Pin
    };

private:
    // MPC definition
    Type type_;
    ID   node_1_id_;
    ID   node_2_id_;

public:
    // Construction
    Mpc(ID node_1_id, ID node_2_id, Type type);

    // Generate the homogeneous linear equations for the selected MPC type.
    // BEAM uses the current nodal positions to evaluate the master-to-slave
    // offset of the rigid-link kinematics.
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const;

    // Additional generalized DOFs required at the independent second node.
    // Existing slave DOFs determine which master translations and rotations
    // must participate in the MPC equations.
    [[nodiscard]] Dofs master_dofs(const SystemDofs& system_dof_mask) const;

    // Node identifiers follow Abaqus MPC ordering: node 1 is dependent and
    // node 2 is the independent master.
    [[nodiscard]] ID node_1() const { return node_1_id_; }
    [[nodiscard]] ID node_2() const { return node_2_id_; }
};

} // namespace constraint
} // namespace fem
