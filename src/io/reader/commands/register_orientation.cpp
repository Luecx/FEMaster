/**
 * @file register_orientation.cpp
 * @brief Registers object-based rectangular and cylindrical coordinate systems.
 *
 * The parser constructs the selected concrete coordinate-system type including
 * its intrinsic name and passes the finished object to
 * `Model::add_coordinate_system()`. This removes the templated coordinate-system
 * factory and duplicate name argument from the model API. Native TYPE-based
 * definitions accept DEFINITION as an unused legacy parameter; Abaqus point-based
 * definitions without TYPE support only DEFINITION=COORDINATES.
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

/**
 * @brief Registers native vector and Abaqus point-based orientation syntax.
 *
 * TYPE selects native rectangular or cylindrical vector definitions, where
 * DEFINITION remains an unused legacy parameter. Without TYPE, SYSTEM selects
 * the supported rectangular Abaqus point definition, with optional
 * DEFINITION=COORDINATES. TYPE and SYSTEM remain mutually exclusive.
 * Executing a matched data variant adds the constructed system to the model.
 *
 * @param registry Registry receiving the orientation command and data variants.
 * @param model Model receiving the named coordinate-system definitions.
 */
void register_orientation(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ORIENTATION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define an Abaqus point-based orientation by default, or explicit FEMaster vector axes.");

        auto name = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional().allowed({"RECTANGULAR", "CYLINDRICAL"})
                .key("SYSTEM").optional().allowed({"RECTANGULAR"})
                .key("DEFINITION").optional().doc("COORDINATES for Abaqus; unused legacy parameter with TYPE")
                .key("NAME").required().doc("Coordinate system identifier")
        );

        command.on_enter([name](const fem::io::dsl::Keys& keys) {
            // TYPE selects native vector semantics, while SYSTEM belongs to the
            // Abaqus point definition. Legacy DEFINITION values do not alter TYPE.
            logging::error(!(keys.has("TYPE") && keys.has("SYSTEM")),
                "ORIENTATION: TYPE and SYSTEM are mutually exclusive");
            logging::error(keys.has("TYPE") || !keys.has("DEFINITION") || keys.raw("DEFINITION") == "COORDINATES",
                "ORIENTATION: without TYPE, DEFINITION must be COORDINATES");

            *name = keys.raw("NAME");
        });

        const auto native_rectangular = fem::io::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"});

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(3))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA").desc("Rectangular system vectors")
                )
                .bind([&model, name](const std::array<fem::Precision, 9>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]},
                        fem::Vec3{values[6], values[7], values[8]}
                    ));
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(2))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 6>().name("DATA").desc("Rectangular system vectors (two vectors)")
                )
                .bind([&model, name](const std::array<fem::Precision, 6>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]}
                    ));
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(native_rectangular)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 3>().name("DATA").desc("Rectangular system vector")
                )
                .bind([&model, name](const std::array<fem::Precision, 3>& values) {
                    model.add_coordinate_system(std::make_shared<cos::RectangularSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]}
                    ));
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"CYLINDRICAL"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(3))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA").desc("Cylindrical system vectors")
                )
                .bind([&model, name](const std::array<fem::Precision, 9>& values) {
                    model.add_coordinate_system(std::make_shared<cos::CylindricalSystem>(
                        *name,
                        fem::Vec3{values[0], values[1], values[2]},
                        fem::Vec3{values[3], values[4], values[5]},
                        fem::Vec3{values[6], values[7], values[8]}
                    ));
                })
            )
        );

        // Abaqus points a, b and optional origin c are the default without TYPE.
        const auto abaqus_rectangular = fem::io::dsl::Condition::any_of({
            fem::io::dsl::Condition::key_equals("SYSTEM", {"RECTANGULAR"}),
            fem::io::dsl::Condition::negate(fem::io::dsl::Condition::key_present("TYPE"))
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(abaqus_rectangular)
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
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
