/**
 * @file register_transform.inl
 * @brief Registers supported Abaqus *TRANSFORM nodal coordinate systems.
 *
 * `*TRANSFORM` associates a rectangular or cylindrical coordinate system with
 * all nodes of one Abaqus node set. The registration creates the corresponding
 * FEMaster coordinate system and stores its name for each affected node so loads
 * and supports can evaluate their components in the transformed basis.
 *
 * Rectangular (`TYPE=R`) and cylindrical (`TYPE=C`) transformations are
 * supported.
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
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../cos/cylindrical_system.h"
#include "../../../cos/rectangular_system.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers rectangular and cylindrical Abaqus nodal transformations.
 *
 * `TYPE=R` uses points `a` and `b` to define local direction 1 and the local
 * 1-2 plane. `TYPE=C` uses `a` and `b` as points on the cylindrical axis; the
 * radial and tangential directions are evaluated by `CylindricalSystem` at each
 * target node. A node may not receive conflicting transform definitions.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing node-to-transform associations.
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

                    logging::error(
                        model._data->node_sets.has(*nset),
                        "TRANSFORM references unknown node set '", *nset, "'"
                    );

                    const fem::Vec3 a{data[0], data[1], data[2]};
                    const fem::Vec3 b{data[3], data[4], data[5]};
                    const fem::Precision eps = std::numeric_limits<fem::Precision>::epsilon();
                    const std::string orientation = "__ABQ_TRANSFORM_" + *nset;

                    logging::error(
                        !model._data->coordinate_systems.has(orientation),
                        "TRANSFORM for node set '", *nset, "' is defined more than once"
                    );

                    if (*type == "R") {
                        const fem::Precision norm_a = a.norm();
                        const fem::Precision norm_b = b.norm();
                        const fem::Precision cross  = a.cross(b).norm();
                        logging::error(
                            norm_a > eps && norm_b > eps && cross > eps * norm_a * norm_b,
                            "TRANSFORM TYPE=R requires nonzero, non-collinear points a and b"
                        );

                        model.add_coordinate_system<cos::RectangularSystem>(orientation, a, b);
                    } else {
                        const fem::Vec3 axis = b - a;
                        logging::error(
                            axis.norm() > eps,
                            "TRANSFORM TYPE=C requires distinct axis points a and b"
                        );

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
                        logging::error(
                            inserted || it->second == orientation,
                            "Node ", node_id,
                            " belongs to multiple incompatible TRANSFORM definitions"
                        );
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
