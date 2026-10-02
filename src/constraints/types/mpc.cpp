/**
 * @file mpc.cpp
 * @brief Implements simple two-node multi-point constraint equations.
 *
 * The implementation converts BEAM, TIE and PIN definitions into the symbolic
 * equations consumed by FEMaster's common constraint system. BEAM follows the
 * same small-rotation rigid-link kinematics as a one-slave kinematic coupling,
 * while TIE and PIN directly equate selected nodal generalized DOFs.
 *
 * @see src/constraints/types/mpc.h
 * @see src/constraints/types/coupling.cpp
 * @see src/constraints/transformer/constraint_system.h
 *
 * @author Finn Eggers
 * @date 01.10.2026
 */

#include "mpc.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

namespace fem {
namespace constraint {

/**
 * Constructs one two-node MPC definition.
 *
 * The first node is retained as the dependent node and the second node as the
 * independent master. No equations are generated until the active structural
 * DOF map is available during model assembly.
 *
 * @param node_1_id Dependent node identifier.
 * @param node_2_id Independent master node identifier.
 * @param type Supported MPC relation.
 */
Mpc::Mpc(ID node_1_id, ID node_2_id, Type type)
    : type_      (type),
      node_1_id_(node_1_id),
      node_2_id_(node_2_id) {}

/**
 * Generates the homogeneous linear equations represented by this MPC.
 *
 * TIE and PIN directly equate generalized nodal DOFs that are active at both
 * nodes. BEAM treats the second node as a rigid-link master. For the offset
 *
 *     r = x_slave - x_master,
 *
 * the translational equations enforce
 *
 *     u_slave = u_master + theta_master x r,
 *
 * and available rotational DOFs satisfy
 *
 *     theta_slave = theta_master.
 *
 * This is the same small-rotation kinematics used by a kinematic coupling with
 * one slave node. The resulting equations remain ordinary rows of C u = 0.
 *
 * @param system_nodal_dofs Active global DOF identifiers.
 * @param model_data Model data containing the current nodal positions.
 * @return Symbolic homogeneous constraint equations.
 */
Equations Mpc::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const {
    Equations equations{};

    // PIN and TIE directly equate the selected active generalized DOFs
    if (type_ == Type::Pin || type_ == Type::Tie) {
        const Dim n_dofs = type_ == Type::Pin ? Dim(3) : Dim(6);

        for (Dim dof = 0; dof < n_dofs; ++dof) {
            const auto dof_id_n1 = system_nodal_dofs(node_1_id_, dof);
            const auto dof_id_n2 = system_nodal_dofs(node_2_id_, dof);
            if (dof_id_n1 < 0 || dof_id_n2 < 0) {
                continue;
            }

            equations.emplace_back(Equation{
                EquationEntry{node_1_id_, dof, Precision(1)},
                EquationEntry{node_2_id_, dof, Precision(-1)}
            });
        }

        return equations;
    }

    // BEAM requires the geometric offset between dependent and master nodes
    logging::error(model_data.positions != nullptr,
        "MPC BEAM: POSITION field is not initialized");

    const auto& node_coords = *model_data.positions;
    const Vec3 r = node_coords.row_vec3(static_cast<Index>(node_1_id_)) - node_coords.row_vec3(static_cast<Index>(node_2_id_));

    const Precision dx = r(0);
    const Precision dy = r(1);
    const Precision dz = r(2);

    // Couple slave translations to the master translation and rotation using
    // u_slave = u_master + theta_master x r
    for (Dim dof = 0; dof < 3; ++dof) {
        if (system_nodal_dofs(node_1_id_, dof) < 0) {
            continue;
        }

        const Dim rotation_1 = static_cast<Dim>((dof + 1) % 3 + 3);
        const Dim rotation_2 = static_cast<Dim>((dof + 2) % 3 + 3);

        Precision coefficient_1 = Precision(0);
        Precision coefficient_2 = Precision(0);

        if (dof == 0) {
            coefficient_1 =  dz;
            coefficient_2 = -dy;
        } else if (dof == 1) {
            coefficient_1 =  dx;
            coefficient_2 = -dz;
        } else {
            coefficient_1 =  dy;
            coefficient_2 = -dx;
        }

        equations.emplace_back(Equation{
            EquationEntry{node_2_id_, dof,        Precision(1)},
            EquationEntry{node_2_id_, rotation_1, coefficient_1},
            EquationEntry{node_2_id_, rotation_2, coefficient_2},
            EquationEntry{node_1_id_, dof,        Precision(-1)}
        });
    }

    // Couple rotational DOFs only when both nodes expose the corresponding
    // generalized component
    for (Dim dof = 3; dof < 6; ++dof) {
        const auto dof_id_n1 = system_nodal_dofs(node_1_id_, dof);
        const auto dof_id_n2 = system_nodal_dofs(node_2_id_, dof);
        if (dof_id_n1 < 0 || dof_id_n2 < 0) {
            continue;
        }

        equations.emplace_back(Equation{
            EquationEntry{node_1_id_, dof, Precision(1)},
            EquationEntry{node_2_id_, dof, Precision(-1)}
        });
    }

    return equations;
}

/**
 * Determines additional active DOFs required at the independent master node.
 *
 * PIN requires only translations already present at the dependent node. TIE
 * mirrors every available dependent generalized DOF. BEAM additionally requires
 * the master rotations that contribute to each active slave translation through
 * the rigid-link cross product.
 *
 * @param system_dof_mask Current structural DOF activation mask.
 * @return Master-node DOFs that must be activated for this MPC.
 */
Dofs Mpc::master_dofs(const SystemDofs& system_dof_mask) const {
    Dofs required{false, false, false, false, false, false};

    // PIN and TIE require a matching master component for every selected slave
    // component that already participates in the structural system
    const Dim n_direct_dofs = type_ == Type::Pin ? Dim(3) : Dim(6);
    for (Dim dof = 0; dof < n_direct_dofs; ++dof) {
        required(dof) = system_dof_mask(node_1_id_, dof);
    }

    if (type_ != Type::Beam) {
        return required;
    }

    // A BEAM master rotation is required whenever it is directly coupled at the
    // slave or contributes to an active slave translation through theta x r
    required(3) = required(3) || system_dof_mask(node_1_id_, 1) || system_dof_mask(node_1_id_, 2);
    required(4) = required(4) || system_dof_mask(node_1_id_, 0) || system_dof_mask(node_1_id_, 2);
    required(5) = required(5) || system_dof_mask(node_1_id_, 0) || system_dof_mask(node_1_id_, 1);

    return required;
}

} // namespace constraint
} // namespace fem
