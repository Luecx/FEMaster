/**
 * @file rigid_body.cpp
 * @brief Region resolution and rigid kinematics for Abaqus-style rigid bodies.
 */

#include "rigid_body.h"

#include "coupling.h"
#include "../../core/logging.h"
#include "../../model/model_data.h"
#include "../../model/element/element_structural.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace fem::constraint {

void RigidBody::prepare(model::ModelData& data) {
    logging::error(data.compiled && data.positions != nullptr,
        "RIGIDBODY: model must be compiled before preparing rigid bodies");
    logging::error(reference_node >= 0 && reference_node < data.positions->rows,
        "RIGIDBODY: invalid reference node ", reference_node);
    logging::error(element_region || pin_region || tie_region,
        "RIGIDBODY: ELSET, PIN NSET or TIE NSET is required");

    // For mixed ELSETs, rotational elements imply TIE unless an explicit
    // PIN/TIE selection overrides the automatically inferred classification.
    std::map<ID, bool> is_tie;
    std::vector<ID> elements;
    if (element_region) {
        logging::error(element_region->size() > 0,
            "RIGIDBODY: empty ELSET ", element_region->name);
        for (ID id : *element_region) {
            logging::error(id >= 0 && static_cast<std::size_t>(id) < data.elements.size()
                && data.elements[static_cast<std::size_t>(id)] != nullptr,
                "RIGIDBODY: invalid element ", id);
            logging::error(!data.is_rigid_element(id),
                "RIGIDBODY: element in multiple rigid bodies: ", id);

            const auto& element = data.elements[static_cast<std::size_t>(id)];
            logging::error(element->as<model::StructuralElement>() != nullptr,
                "RIGIDBODY: element is not structural: ", id);
            const ElDofs dofs = element->dofs();
            const bool rotates = dofs(0, 3) || dofs(0, 4) || dofs(0, 5);
            elements.push_back(id);
            for (ID node : *element) {
                auto [it, inserted] = is_tie.emplace(node, rotates);
                if (!inserted) it->second = it->second || rotates;
            }
        }
        std::sort(elements.begin(), elements.end());
        elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
    }

    std::set<ID> explicit_pin, explicit_tie;
    auto copy_explicit = [&](const model::NodeRegion::Ptr& region, std::set<ID>& target) {
        if (!region) return;
        logging::error(region->size() > 0,
            "RIGIDBODY: empty node set ", region->name);
        for (ID node : *region) {
            logging::error(node >= 0 && node < data.positions->rows,
                "RIGIDBODY: invalid node ", node);
            target.insert(node);
        }
    };
    copy_explicit(pin_region, explicit_pin);
    copy_explicit(tie_region, explicit_tie);

    for (ID node : explicit_pin) {
        logging::error(explicit_tie.count(node) == 0,
            "RIGIDBODY: node in both PIN NSET and TIE NSET: ", node);
        is_tie[node] = false;
    }
    for (ID node : explicit_tie) is_tie[node] = true;

    // Prevent overlapping rigid bodies and a reference node that is the slave
    // of another body. A reference node never constrains itself.
    for (const auto& previous : data.rigid_bodies) {
        logging::error(previous.reference_node != reference_node,
            "RIGIDBODY: reference node reused by another rigid body");
        for (const auto& [node, tie] : is_tie) {
            (void)tie;
            if (node == reference_node) continue;
            logging::error(node != previous.reference_node,
                "RIGIDBODY: node is the reference of another rigid body: ", node);
            const bool shared_pin = previous.pin_nodes()
                && previous.pin_nodes()->contains(node);
            const bool shared_tie = previous.tie_nodes()
                && previous.tie_nodes()->contains(node);
            logging::error(!shared_pin && !shared_tie,
                "RIGIDBODY: slave node belongs to multiple rigid bodies: ", node);
        }
        const bool ref_is_pin = previous.pin_nodes()
            && previous.pin_nodes()->contains(reference_node);
        const bool ref_is_tie = previous.tie_nodes()
            && previous.tie_nodes()->contains(reference_node);
        logging::error(!ref_is_pin && !ref_is_tie,
            "RIGIDBODY: reference node belongs to another rigid body");
    }

    auto pin = std::make_shared<model::NodeRegion>(
        "RIGIDBODY_PIN_" + std::to_string(reference_node));
    auto tie = std::make_shared<model::NodeRegion>(
        "RIGIDBODY_TIE_" + std::to_string(reference_node));
    for (const auto& [node, rotates] : is_tie) {
        if (node == reference_node) continue;
        if (rotates) tie->add(node);
        else         pin->add(node);
    }

    if (position == Position::CenterOfMass) {
        Precision total_mass = Precision(0);
        Precision reference_mass = Precision(0);
        Vec3 first_moment = Vec3::Zero();

        auto accumulate = [&](const model::ElementPtr& element) {
            auto* structural = element->as<model::StructuralElement>();
            if (!structural) return;
            const Precision m = structural->integrate_scalar_field(
                true, [](const Vec3&) { return Precision(1); });
            logging::error(std::isfinite(m) && m >= Precision(0),
                "RIGIDBODY: invalid element mass");
            if (m == Precision(0)) return;
            const Vec3 moment = structural->integrate_vector_field(
                true, [](const Vec3& x) { return x; });
            logging::error(moment.allFinite(), "RIGIDBODY: invalid first mass moment");
            total_mass += m;
            first_moment += moment;
        };
        for (ID id : elements) {
            accumulate(data.elements[static_cast<std::size_t>(id)]);
        }

        // Both regular MASS elements and auxiliary POINTMASS objects can
        // contribute even without being selected by ELSET. The reference
        // node's own point mass does not change its initial location.
        auto add_point = [&](const model::ElementPtr& element) {
            if (!element || element->n_nodes() != 1) return;
            const ID node = element->nodes()[0];
            if (node == reference_node) {
                auto* structural = element->as<model::StructuralElement>();
                if (structural) reference_mass += structural->integrate_scalar_field(
                    true, [](const Vec3&) { return Precision(1); });
            } else if (is_tie.count(node) != 0) {
                accumulate(element);
            }
        };
        for (const auto& element : data.elements) {
            if (!element || element->n_nodes() != 1
                || std::binary_search(elements.begin(), elements.end(), element->elem_id)) continue;
            add_point(element);
        }
        for (const auto& point : data.point_elements) add_point(point);

        logging::error(std::isfinite(total_mass) && total_mass >= Precision(0)
            && std::isfinite(reference_mass) && reference_mass >= Precision(0)
            && first_moment.allFinite(),
            "RIGIDBODY: invalid total mass for CENTER OF MASS");
        if (total_mass == Precision(0)) {
            logging::error(reference_mass > Precision(0),
                "RIGIDBODY: CENTER OF MASS requires positive mass");
        } else {
            const Vec3 center = first_moment / total_mass;
            for (Index d = 0; d < 3; ++d) {
                (*data.positions)(reference_node, d) = center[d];
                if (data.positions_reference && data.positions_reference != data.positions) {
                    (*data.positions_reference)(reference_node, d) = center[d];
                }
            }
        }
    }

    // Commit only once region resolution and mass validation have succeeded.
    elements_ = std::move(elements);
    resolved_pin_ = std::move(pin);
    resolved_tie_ = std::move(tie);
    if (data.rigid_element_mask.size() < data.elements.size()) {
        data.rigid_element_mask.resize(data.elements.size(), 0);
    }
    for (ID id : elements_) {
        data.rigid_element_mask[static_cast<std::size_t>(id)] = 1;
    }
}

Equations RigidBody::get_equations(SystemDofIds& indices, model::ModelData& data) const {
    Equations all;
    auto append = [&](const model::NodeRegion::Ptr& nodes, const Dofs& dofs) {
        if (!nodes || nodes->size() == 0) return;
        Coupling coupling{reference_node, nodes, dofs, CouplingType::KINEMATIC};
        auto part = coupling.get_equations(indices, data);
        all.insert(all.end(), std::make_move_iterator(part.begin()),
                             std::make_move_iterator(part.end()));
    };
    append(resolved_pin_, Dofs{true, true, true, false, false, false});
    append(resolved_tie_, Dofs{true, true, true, true, true, true});
    return all;
}

Dofs RigidBody::master_dofs(SystemDofs& mask, model::ModelData& data) const {
    Dofs necessary{false, false, false, false, false, false};
    auto include = [&](const model::NodeRegion::Ptr& nodes, const Dofs& dofs) {
        if (!nodes || nodes->size() == 0) return;
        Coupling coupling{reference_node, nodes, dofs, CouplingType::KINEMATIC};
        const Dofs active = coupling.master_dofs(mask, data);
        for (Index d = 0; d < 6; ++d) necessary(d) = necessary(d) || active(d);
    };
    include(resolved_pin_, Dofs{true, true, true, false, false, false});
    include(resolved_tie_, Dofs{true, true, true, true, true, true});
    return necessary;
}

} // namespace fem::constraint
