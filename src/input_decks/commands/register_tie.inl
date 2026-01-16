// register_tie.inl — registers *TIE

#include <stdexcept>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_tie(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("TIE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Bind a slave set (node or surface set) to a master surface/line set (tie constraint).");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MASTER")  .required(    ).doc("Master surface set or line set")
                .key("SLAVE")   .required(    ).doc("Slave node set or slave surface set")
                .key("ADJUST")  .optional("NO").doc("Determines projection of slave to master").allowed({"NO", "YES"})
                .key("DISTANCE").required(    ).doc("Search distance")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            const std::string& master = keys.raw("MASTER");
            const std::string& slave  = keys.raw("SLAVE");
            const bool adjust         = keys.get<bool>("ADJUST");
            const fem::Precision dist = keys.get<fem::Precision>("DISTANCE");
            model.add_tie(master, slave, dist, adjust);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
