// register_beam_section.inl — DSL registration for *BEAMSECTION

#include <array>
#include <memory>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_beam_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("BEAMSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a beam section with orientation to an element set.");

        // Persistent per-command state (lives beyond this function):
        auto material = std::make_shared<std::string>();
        auto elset    = std::make_shared<std::string>();
        auto profile  = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required().doc("Material name")
                .key("ELSET").required().doc("Target element set")
                .key("PROFILE").required().doc("Beam profile identifier")
        );

        // Capture shared_ptrs BY VALUE so they remain valid when executed later
        command.on_enter([material, elset, profile](const fem::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
            *profile  = keys.raw("PROFILE");
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<fem::Precision, 3>().name("N1").desc("Section orientation vector")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, material, elset, profile](const std::array<fem::Precision, 3>& n1_data) {
                    fem::Vec3 n1; n1 << n1_data[0], n1_data[1], n1_data[2];
                    model.beam_section(*elset, *material, *profile, n1);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
