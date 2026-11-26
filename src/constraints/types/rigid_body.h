/**
 * @file rigid_body.h
 * @brief Declares the Abaqus-style rigid body constraint.
 *
 * A rigid body ties the nodes of an element set and of optional tie and pin
 * node sets to one reference node. Tie nodes follow translations and rotations
 * of the reference node, pin nodes only its translations. Nodes of the element
 * set are tie nodes unless they appear in the pin node set; rotations are only
 * constrained where the slave node carries active rotational DOFs.
 *
 * The equations are generated through kinematic couplings so the rigid-body
 * kinematics exist in one place only.
 *
 * @see Coupling
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#pragma once

#include "coupling.h"

#include "../../core/printable.h"
#include "../../data/region.h"
#include "../../model/system_dofs.h"
#include "equation.h"

#include <string>
#include <vector>

namespace fem {
namespace constraint {

class RigidBody : public fem::Printable {
public:
    Index                     ref_node;
    model::NodeRegion::Ptr    ref_region = nullptr;
    model::ElementRegion::Ptr elements   = nullptr;
    model::NodeRegion::Ptr    tie_nodes  = nullptr;
    model::NodeRegion::Ptr    pin_nodes  = nullptr;

public:
    // At least one of the element, tie or pin sets must be given. A node must
    // not appear in both the tie and the pin node set.
    RigidBody(Index                     ref_node,
              model::ElementRegion::Ptr elements,
              model::NodeRegion::Ptr    tie_nodes,
              model::NodeRegion::Ptr    pin_nodes);

    Equations   get_equations(model::SystemDofs& system_nodal_dofs, model::ModelData& model_data);
    Dofs        master_dofs  (const model::SystemDofs& system_dofs, model::ModelData& model_data);
    std::string str() const override;

private:
    // Kinematic couplings for the tie nodes (all six DOFs) and pin nodes
    // (translations only) of this rigid body
    std::vector<Coupling> couplings(model::ModelData& model_data) const;
};

} // namespace constraint
} // namespace fem
