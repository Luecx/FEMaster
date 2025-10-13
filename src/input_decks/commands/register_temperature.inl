// register_temperature.inl — registers *TEMPERATURE

#include <stdexcept>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_temperature(fem::dsl::Registry& registry, model::Model& model) {
    std::string field_name;

    registry.command("TEMPERATURE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign nodal temperatures (supports sets or individual nodes).");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Temperature field name")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            field_name = keys.raw("NAME");
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or id")
                    .one<fem::Precision>().name("VALUE").desc("Temperature value")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&](const std::string& target, fem::Precision value) {
                    if (model._data->node_sets.has(target)) {
                        for (auto id : *model._data->node_sets.get(target)) {
                            model.set_field_temperature(field_name, id, value);
                        }
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.set_field_temperature(field_name, id, value);
                    } catch (const std::exception&) {
                        throw std::runtime_error("TEMPERATURE target '" + target + "' is neither a node set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
