/**
 * @file register_orientation.inl
 * @brief Registers coordinate-defined rectangular Abaqus *ORIENTATION data.
 *
 * The initial Abaqus orientation reader supports the default
 * `DEFINITION=COORDINATES` and `SYSTEM=RECTANGULAR` form. Abaqus points `a`, `b`
 * and optional origin `c` are converted to the two direction vectors required
 * by FEMaster's `RectangularSystem`.
 *
 * FEMaster coordinate systems currently represent orientation only and do not
 * store a translational origin. This is sufficient for material and section
 * orientations because only the basis directions are consumed. The Abaqus
 * origin therefore contributes to the direction vectors `a-c` and `b-c` but is
 * not retained afterwards.
 *
 * Node-defined, offset-to-node, cylindrical, spherical, Z-rectangular, user and
 * additional-rotation orientation forms remain unsupported.
 *
 * @see cos::RectangularSystem
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../cos/rectangular_system.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the basic rectangular Abaqus orientation definition.
 *
 * The first six data values define points `a` and `b`. The optional final three
 * values define point `c`; omitted values place `c` at the global origin. The
 * local 1-direction is constructed from `a-c`, while `b-c` defines the local
 * 1-2 plane. `RectangularSystem` performs the final orthogonalization.
 *
 * Degenerate input is rejected when `a` coincides with `c` or when the two
 * supplied directions are collinear. A second Abaqus data line for an additional
 * rotation is intentionally not accepted by this initial implementation.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the coordinate system.
 */
inline void register_orientation(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ORIENTATION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define a rectangular Abaqus material orientation from points a, b and optional c.");

        // Preserve the orientation name until the following data row constructs
        // the coordinate system
        auto name = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Orientation identifier")
                .key("DEFINITION")
                    .optional("COORDINATES")
                    .allowed({"COORDINATES"})
                    .doc("Orientation definition method")
                .key("SYSTEM")
                    .optional("RECTANGULAR")
                    .allowed({"RECTANGULAR"})
                    .doc("Local coordinate-system type")
        );

        command.on_enter([name](const fem::io::dsl::Keys& keys) {
            *name = keys.raw("NAME");
        });

        // Read the coordinate-defined Abaqus points on exactly one data line.
        // Point c defaults to the global origin when its coordinates are omitted.
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("A1").desc("Point a, global x")
                    .one<fem::Precision>().name("A2").desc("Point a, global y")
                    .one<fem::Precision>().name("A3").desc("Point a, global z")
                    .one<fem::Precision>().name("B1").desc("Point b, global x")
                    .one<fem::Precision>().name("B2").desc("Point b, global y")
                    .one<fem::Precision>().name("B3").desc("Point b, global z")
                    .one<fem::Precision>().name("C1").desc("Point c, global x")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("C2").desc("Point c, global y")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("C3").desc("Point c, global z")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, name](fem::Precision a1,
                                     fem::Precision a2,
                                     fem::Precision a3,
                                     fem::Precision b1,
                                     fem::Precision b2,
                                     fem::Precision b3,
                                     fem::Precision c1,
                                     fem::Precision c2,
                                     fem::Precision c3) {
                    const fem::Vec3 a{a1, a2, a3};
                    const fem::Vec3 b{b1, b2, b3};
                    const fem::Vec3 c{c1, c2, c3};

                    // Convert Abaqus point coordinates to the direction vectors
                    // represented by FEMaster's origin-free coordinate system
                    const fem::Vec3 axis_1   = a - c;
                    const fem::Vec3 in_plane = b - c;

                    const fem::Precision norm_1     = axis_1.norm();
                    const fem::Precision norm_plane = in_plane.norm();
                    const fem::Precision cross_norm = axis_1.cross(in_plane).norm();
                    const fem::Precision tolerance  = std::numeric_limits<fem::Precision>::epsilon();

                    if (norm_1 <= tolerance || norm_plane <= tolerance ||
                        cross_norm <= tolerance * norm_1 * norm_plane) {
                        throw std::runtime_error("ORIENTATION requires distinct, non-collinear points a, b and c");
                    }

                    // RectangularSystem orthogonalizes the in-plane direction
                    // against axis 1 and completes the right-handed basis
                    model.add_coordinate_system<cos::RectangularSystem>(*name, axis_1, in_plane);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
