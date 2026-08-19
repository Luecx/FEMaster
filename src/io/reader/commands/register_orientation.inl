/**
 * @file register_orientation.inl
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

#include <array>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../cos/cylindrical_system.h"
#include "../../../cos/rectangular_system.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_orientation(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ORIENTATION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define a rectangular or cylindrical coordinate system.");

        auto name = std::make_shared<std::string>();
        auto type = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"RECTANGULAR", "CYLINDRICAL"})
                .key("DEFINITION").optional().doc("Unused legacy parameter")
                .key("NAME").required().doc("Coordinate system identifier")
        );

        command.on_enter([name, type](const fem::io::dsl::Keys& keys) {
            *type = keys.raw("TYPE");
            *name = keys.raw("NAME");
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
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
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
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
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
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
    });
}

} // namespace fem::io::reader::commands