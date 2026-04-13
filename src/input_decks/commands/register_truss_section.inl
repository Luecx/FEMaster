// register_truss_section.inl — DSL registration for *TRUSSSECTION

#include <memory>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_truss_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("TRUSSSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a truss section area to an element set.");

        auto material = std::make_shared<std::string>();
        auto elset    = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required().doc("Material name")
                .key("ELSET").required().doc("Target truss element set")
        );

        command.on_enter([material, elset](const fem::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("AREA").desc("Truss cross-sectional area")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, material, elset](fem::Precision area) {
                    model.truss_section(*elset, *material, area);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
