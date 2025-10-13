// register_orientation.inl — registers *ORIENTATION (coordinate systems)

#include <array>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../cos/cylindrical_system.h"
#include "../../cos/rectangular_system.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_orientation(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("ORIENTATION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define coordinate systems (rectangular or cylindrical).");

        // Persistent per-command state (lives beyond this function)
        auto name = std::make_shared<std::string>();
        auto type = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"RECTANGULAR", "CYLINDRICAL"})
                .key("DEFINITION").optional().doc("Unused legacy parameter")
                .key("NAME").required().doc("Coordinate system identifier")
        );

        // Capture shared_ptrs BY VALUE so they're safe later
        command.on_enter([name, type](const fem::dsl::Keys& keys) {
            *type = keys.raw("TYPE");
            *name = keys.raw("NAME");
        });

        // RECTANGULAR variants -------------------------------------------------
        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA").desc("Rectangular system vectors")
                )
                .bind([&model, name](const std::array<fem::Precision, 9>& vals) {
                    model.add_coordinate_system<cos::RectangularSystem>(*name,
                        fem::Vec3{vals[0], vals[1], vals[2]},
                        fem::Vec3{vals[3], vals[4], vals[5]},
                        fem::Vec3{vals[6], vals[7], vals[8]});
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 6>().name("DATA").desc("Rectangular system vectors (two vectors)")
                )
                .bind([&model, name](const std::array<fem::Precision, 6>& vals) {
                    model.add_coordinate_system<cos::RectangularSystem>(*name,
                        fem::Vec3{vals[0], vals[1], vals[2]},
                        fem::Vec3{vals[3], vals[4], vals[5]});
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"RECTANGULAR"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 3>().name("DATA").desc("Rectangular system vector")
                )
                .bind([&model, name](const std::array<fem::Precision, 3>& vals) {
                    model.add_coordinate_system<cos::RectangularSystem>(*name,
                        fem::Vec3{vals[0], vals[1], vals[2]});
                })
            )
        );

        // CYLINDRICAL ----------------------------------------------------------
        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"CYLINDRICAL"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 9>().name("DATA").desc("Cylindrical system vectors")
                )
                .bind([&model, name](const std::array<fem::Precision, 9>& vals) {
                    model.add_coordinate_system<cos::CylindricalSystem>(*name,
                        fem::Vec3{vals[0], vals[1], vals[2]},
                        fem::Vec3{vals[3], vals[4], vals[5]},
                        fem::Vec3{vals[6], vals[7], vals[8]});
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
