// register_loadcase_topoorient.inl — registers TOPOORIENT for LINEARSTATICTOPO loadcases

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_topoorient(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPOORIENT", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Select the orientation field for LINEARSTATICTOPO loadcases.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("Orientation field name (ELEMENT, 3 components)")
                .alternative("NAME")
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* lc = parser.active_loadcase_as<loadcase::LinearStaticTopo>();
            if (!lc) {
                throw std::runtime_error("TOPOORIENT only valid for LINEARSTATICTOPO loadcases");
            }

            const std::string field_name = keys.raw("FIELD");
            auto field = parser.model()._data->get_field(field_name);
            if (!field) {
                throw std::runtime_error("TOPOORIENT field '" + field_name + "' does not exist");
            }
            if (field->domain != model::FieldDomain::ELEMENT || field->components != 3) {
                throw std::runtime_error("TOPOORIENT field '" + field_name + "' must be ELEMENT domain with 3 components");
            }
            lc->orientation = field;
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(0).max(0))
                .pattern(fem::dsl::Pattern::make())
                .bind([]() {})
            )
        );
    });
}

} // namespace fem::input_decks::commands
