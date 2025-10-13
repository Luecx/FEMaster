// ===== FILE: ./include/fem/input_decks/commands/register_profile.inl =====
#include <array>
#include <memory>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_profile(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("PROFILE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define a beam profile (A, Iy, Iz, Jt).");

        // Persistent state for this command's handlers:
        auto profile_name = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("PROFILE")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the profile")
        );

        // Capture by value (shared_ptr) so it stays valid
        command.on_enter([profile_name](const fem::dsl::Keys& keys) {
            *profile_name = keys.raw("PROFILE");  // spec canonicalizes NAME→PROFILE
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<fem::Precision, 4>().name("DATA").desc("Area, Iy, Iz, Jt")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, profile_name](const std::array<fem::Precision, 4>& data) {
                    // use *profile_name safely
                    model._data->profiles.activate(*profile_name, data[0], data[1], data[2], data[3]);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
