// register_contact.inl - registers *CONTACT

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_contact(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CONTACT", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define frictionless node-to-surface penalty contact. MASTER must be a surface set; "
            "SLAVE may be a node set or surface set. Contact contributes only in NONLINEARSTATIC."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MASTER")   .required()     .doc("Master surface set")
                .key("SLAVE")    .required()     .doc("Slave node set or slave surface set")
                .key("DISTANCE") .required()     .doc("Search distance for candidate master surfaces")
                .key("PENALTY")  .required()     .doc("Normal penalty stiffness")
                .key("CLEARANCE").optional("0")  .doc("Allowed normal clearance before contact activates")
                .key("FLIP")     .optional("NO") .doc("Flip master surface normals").allowed({"NO", "YES"})
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            const std::string& master = keys.raw("MASTER");
            const std::string& slave  = keys.raw("SLAVE");
            const fem::Precision dist = keys.get<fem::Precision>("DISTANCE");
            const fem::Precision k    = keys.get<fem::Precision>("PENALTY");
            const fem::Precision c    = keys.get<fem::Precision>("CLEARANCE");
            const bool flip           = keys.get<bool>("FLIP");

            model.add_contact(master, slave, dist, k, c, flip);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
