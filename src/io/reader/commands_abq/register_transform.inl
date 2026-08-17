/**
 * @file register_transform.inl
 * @brief Registers supported Abaqus *TRANSFORM nodal coordinate systems.
 *
 * Abaqus associates a transformation with the nodes in one node set. FEMaster
 * intentionally does not transform nodal degrees of freedom globally; loads and
 * supports already carry an optional coordinate system and evaluate its basis at
 * each target node. This registration therefore creates one internal FEMaster
 * coordinate system and records its name for every node in the Abaqus NSET.
 *
 * Rectangular (`TYPE=R`) and cylindrical (`TYPE=C`) transformations are
 * supported. Spherical transformations remain unsupported because FEMaster has
 * no corresponding spatial coordinate-system implementation.
 *
 * @see ParserAbqState
 * @see cos::RectangularSystem
 * @see cos::CylindricalSystem
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../cos/cylindrical_system.h"
#include "../../../cos/rectangular_system.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus rectangular and cylindrical nodal transformations.
 *
 * `TYPE=R` uses point `a` to define local direction 1 and point `b` to define
 * the local 1-2 plane. `TYPE=C` uses points `a` and `b` as two points on the
 * cylindrical axis. For the cylindrical case an arbitrary orthogonal reference
 * radial direction is constructed; the physical radial/tangential basis is then
 * evaluated independently at each node by `CylindricalSystem`.
 *
 * A node may belong to several sets, but it may not receive conflicting Abaqus
 * transformations. Repeating the identical association is harmless.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining node-to-transform associations.
 */
inline void register_transform(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("TRANSFORM", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a rectangular or cylindrical Abaqus transform to a node set.");

        auto nset = std::make_shared<std::string>();
        auto type = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET")
                    .required()
                    .doc("Node set receiving the transformation")
                .key("TYPE")
                    .optional("R")
                    .allowed({"R", "C"})
                    .doc("R = rectangular, C = cylindrical")
        );

        command.on_enter([nset, type](const fem::io::dsl::Keys& keys) {
            *nset = keys.raw("NSET");
            *type = keys.raw("TYPE");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 6>().name("DATA")
                        .desc("Coordinates of Abaqus points a and b")
                )
                .bind([&parser, nset, type](const std::array<fem::Precision, 6>& data) {
                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();

                    if (!model._data->node_sets.has(*nset)) {
                        throw std::runtime_error("TRANSFORM references unknown node set '" + *nset + "'");
                    }

                    const fem::Vec3 a{data[0], data[1], data[2]};
                    const fem::Vec3 b{data[3], data[4], data[5]};
                    const fem::Precision eps = std::numeric_limits<fem::Precision>::epsilon();
                    const std::string orientation = "__ABQ_TRANSFORM_" + *nset;

                    if (model._data->coordinate_systems.has(orientation)) {
                        throw std::runtime_error("TRANSFORM for node set '" + *nset + "' is defined more than once");
                    }

                    if (*type == "R") {
                        const fem::Precision norm_a = a.norm();
                        const fem::Precision norm_b = b.norm();
                        const fem::Precision cross  = a.cross(b).norm();
                        if (norm_a <= eps || norm_b <= eps || cross <= eps * norm_a * norm_b) {
                            throw std::runtime_error("TRANSFORM TYPE=R requires nonzero, non-collinear points a and b");
                        }

                        model.add_coordinate_system<cos::RectangularSystem>(orientation, a, b);
                    } else {
                        const fem::Vec3 axis = b - a;
                        if (axis.norm() <= eps) {
                            throw std::runtime_error("TRANSFORM TYPE=C requires distinct axis points a and b");
                        }

                        const fem::Vec3 axial      = axis.normalized();
                        const fem::Vec3 radial     = axial.unitOrthogonal();
                        const fem::Vec3 tangential = axial.cross(radial).normalized();

                        model.add_coordinate_system<cos::CylindricalSystem>(
                            orientation,
                            a,
                            a + radial,
                            a + tangential
                        );
                    }

                    auto nodes = model._data->node_sets.get(*nset);
                    for (const fem::ID node_id : *nodes) {
                        auto [it, inserted] = state.node_transforms.emplace(node_id, orientation);
                        if (!inserted && it->second != orientation) {
                            throw std::runtime_error(
                                "Node " + std::to_string(node_id) +
                                " belongs to multiple incompatible TRANSFORM definitions"
                            );
                        }
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
