// register_tload.inl — DSL registration for *TLOAD

#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_tload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("TLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Thermal load referencing a temperature field.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .required()
                    .doc("Target load collector that groups the loads")
                .key("TEMPERATUREFIELD")
                    .required()
                    .doc("Name of the temperature field to apply")
                .key("REFERENCETEMPERATURE")
                    .required()
                    .doc("Reference temperature value")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            std::string field            = keys.raw("TEMPERATUREFIELD");
            const fem::Precision ref     = keys.get<fem::Precision>("REFERENCETEMPERATURE");

            model._data->load_cols.activate(collector);
            model.add_tload(field, ref);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
