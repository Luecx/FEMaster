/**
 * @file register_orientation.cpp
 * @brief Registers object-based rectangular and cylindrical coordinate systems.
 *
 * The parser constructs the selected concrete coordinate-system type including
 * its intrinsic name and passes the finished object to
 * `Model::add_coordinate_system()`. This removes the templated coordinate-system
 * factory and duplicate name argument from the model API.
 *
 * @see cos::RectangularSystem
 * @see cos::CylindricalSystem
 * @see model::Model::add_coordinate_system
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <limits>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../core/logging.h"
#include "../../../cos/cylindrical_system.h"
#include "../../../cos/rectangular_system.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

void register_orientation(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ORIENTATION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define an Abaqus point-based orientation by default, or explicit FEMaster vector axes.");

        auto name = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional().allowed({"RECTANGULAR", "CYLINDRICAL"})
                .key("SYSTEM").optional().allowed({"RECTANGULAR"})
                .key("DEFINITION").optional().allowed({"COORDINATES"})
                .key("NAME").required().doc("Coordinate system identifier")
        );

        command.on_enter([name](const fem::io::dsl::Keys& keys) {
            logging::error(!(keys.has("TYPE") && keys.has("SYSTEM")),
                "ORIENTATION: TYPE and SYSTEM are mutually exclusive");
            logging::error(!keys.has("DEFINITION") || !keys.has("TYPE"),
                "ORIENTATION: DEFINITION=COORDINATES is incompatible with TYPE");

            *name = keys.raw("NAME");
        });

        const auto native_rectangular = fem::io::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"});

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .data(
                fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>("DATA", "Rectangular system vectors"),
                [&model, name](const std::array<fem::Precision, 9>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]},
                        fem::Vec3{values[6], values[7], values[8]}
                    ));
                },
                fem::io::dsl::LineRange{}.min(1).max(3)
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .data(
                fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 6>("DATA", "Rectangular system vectors (two vectors)"),
                [&model, name](const std::array<fem::Precision, 6>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]}
                    ));
                },
                fem::io::dsl::LineRange{}.min(1).max(2)
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .data(
                fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 3>("DATA", "Rectangular system vector"),
                [&model, name](const std::array<fem::Precision, 3>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]}
                    ));
                },
                fem::io::dsl::LineRange{}.min(1).max(1)
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"CYLINDRICAL"}))
            .data(
                fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>("DATA", "Cylindrical system vectors"),
                [&model, name](const std::array<fem::Precision, 9>& values) {
                    model.add_coordinate_system(std::make_shared<cos::CylindricalSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]},
                        fem::Vec3{values[6], values[7], values[8]}
                    ));
                },
                fem::io::dsl::LineRange{}.min(1).max(3)
            )
        );

        // Abaqus points a, b and optional origin c are the default without TYPE.
        const auto abaqus_rectangular = fem::io::dsl::Condition::any_of({
            fem::io::dsl::Condition::key_equals("SYSTEM", {"RECTANGULAR"}),
            fem::io::dsl::Condition::negate(fem::io::dsl::Condition::key_present("TYPE"))
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(abaqus_rectangular)
            .data(
                fem::io::dsl::Pattern::make()
                    .one<fem::Precision>("A1", "Point a, global x")
                    .one<fem::Precision>("A2", "Point a, global y")
                    .one<fem::Precision>("A3", "Point a, global z")
                    .one<fem::Precision>("B1", "Point b, global x")
                    .one<fem::Precision>("B2", "Point b, global y")
                    .one<fem::Precision>("B3", "Point b, global z")
                    .one<fem::Precision>("C1", "Point c, global x")
                        .defaults(fem::Precision{0})
                    .one<fem::Precision>("C2", "Point c, global y")
                        .defaults(fem::Precision{0})
                    .one<fem::Precision>("C3", "Point c, global z")
                        .defaults(fem::Precision{0}),
                [&model, name](fem::Precision a1,
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

                    const fem::Vec3 axis_1   = a - c;
                    const fem::Vec3 in_plane = b - c;

                    const fem::Precision norm_1     = axis_1.norm();
                    const fem::Precision norm_plane = in_plane.norm();
                    const fem::Precision cross_norm = axis_1.cross(in_plane).norm();
                    const fem::Precision tolerance  = std::numeric_limits<fem::Precision>::epsilon();

                    logging::error(norm_1 > tolerance
                                && norm_plane > tolerance
                                && cross_norm > tolerance * norm_1 * norm_plane,
                        "ORIENTATION requires distinct, non-collinear points a, b and c");

                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        axis_1,
                        in_plane
                    ));
                },
                fem::io::dsl::LineRange{}.min(1).max(1)
            )
        );
    });
}

} // namespace fem::io::reader::commands
