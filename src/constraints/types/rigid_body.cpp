/**
 * @file rigid_body.cpp
 * @brief Implements the Abaqus-style rigid body constraint.
 *
 * The element set is reduced to its nodes. Together with the tie node set they
 * form the tie nodes; nodes of the pin node set are removed from them and only
 * follow the reference-node translations.
 *
 * @see RigidBody
 * @see Coupling
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#include "rigid_body.h"

#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <sstream>
#include <utility>

namespace fem {
namespace constraint {

/**
 * Constructs a rigid body from its reference node and slave sets.
 *
 * @param ref_node Compiled reference node.
 * @param elements Optional element set whose nodes become tie nodes.
 * @param tie_nodes Optional node set following translations and rotations.
 * @param pin_nodes Optional node set following translations only.
 */
RigidBody::RigidBody(Index                     ref_node,
                     model::ElementRegion::Ptr elements,
                     model::NodeRegion::Ptr    tie_nodes,
                     model::NodeRegion::Ptr    pin_nodes)
    : ref_node (ref_node),
      elements (std::move(elements)),
      tie_nodes(std::move(tie_nodes)),
      pin_nodes(std::move(pin_nodes)) {
    logging::error(this->elements || this->tie_nodes || this->pin_nodes,
        "RigidBody: requires an element set, a tie node set or a pin node set");

    if (this->tie_nodes && this->pin_nodes) {
        for (Index node_id : *this->pin_nodes) {
            logging::error(std::find(this->tie_nodes->begin(), this->tie_nodes->end(), node_id) == this->tie_nodes->end(),
                "RigidBody: node ", node_id, " is in both the tie and the pin node set");
        }
    }
}

/**
 * Builds the kinematic couplings of this rigid body.
 *
 * Tie nodes are the nodes of the element set and the tie node set without the
 * pin nodes. Empty node lists produce no coupling.
 */
std::vector<Coupling> RigidBody::couplings(model::ModelData& model_data) const {
    std::vector<Index> pin_ids;
    if (pin_nodes) {
        pin_ids.assign(pin_nodes->begin(), pin_nodes->end());
    }
    std::sort(pin_ids.begin(), pin_ids.end());

    std::vector<Index> tie_ids;
    if (tie_nodes) {
        tie_ids.assign(tie_nodes->begin(), tie_nodes->end());
    }
    if (elements) {
        for (Index elem_id : *elements) {
            const auto& element = model_data.elements[elem_id];
            logging::error(element != nullptr,
                "RigidBody: element ", elem_id, " of set ", elements->name, " does not exist");
            for (Index node_id : *element) {
                tie_ids.push_back(node_id);
            }
        }
    }
    std::sort(tie_ids.begin(), tie_ids.end());
    tie_ids.erase(std::unique(tie_ids.begin(), tie_ids.end()), tie_ids.end());

    // Pin node types override the tie type of element-set nodes
    auto tie_region = std::make_shared<model::NodeRegion>("RIGIDBODY_TIE");
    for (Index node_id : tie_ids) {
        if (!std::binary_search(pin_ids.begin(), pin_ids.end(), node_id)) {
            tie_region->add(node_id);
        }
    }

    std::vector<Coupling> result;
    if (tie_region->size() > 0) {
        result.emplace_back(ref_node, tie_region,
                            Dofs{true, true, true, true, true, true}, CouplingType::KINEMATIC);
    }
    if (pin_nodes && pin_nodes->size() > 0) {
        result.emplace_back(ref_node, pin_nodes,
                            Dofs{true, true, true, false, false, false}, CouplingType::KINEMATIC);
    }
    return result;
}

Equations RigidBody::get_equations(model::SystemDofs& system_nodal_dofs, model::ModelData& model_data) {
    Equations equations;
    for (Coupling& coupling : couplings(model_data)) {
        Equations coupling_equations = coupling.get_equations(system_nodal_dofs, model_data);
        equations.insert(equations.end(), coupling_equations.begin(), coupling_equations.end());
    }
    return equations;
}

Dofs RigidBody::master_dofs(const model::SystemDofs& system_dofs, model::ModelData& model_data) {
    Dofs dofs{false, false, false, false, false, false};
    for (Coupling& coupling : couplings(model_data)) {
        const Dofs coupling_dofs = coupling.master_dofs(system_dofs, model_data);
        for (Index i = 0; i < 6; ++i) {
            dofs(i) = dofs(i) || coupling_dofs(i);
        }
    }
    return dofs;
}

std::string RigidBody::str() const {
    std::ostringstream s;
    s << "Rigid body: ref=";
    if (ref_region) {
        s << ref_region->name << " (node " << ref_node << ")";
    } else {
        s << "node " << ref_node;
    }
    if (elements)  s << ", ELSET "    << elements->name;
    if (tie_nodes) s << ", TIE NSET " << tie_nodes->name;
    if (pin_nodes) s << ", PIN NSET " << pin_nodes->name;
    return s.str();
}

} // namespace constraint
} // namespace fem
