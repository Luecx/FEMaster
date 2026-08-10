// register_contact.inl - registers *CONTACT

#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_contact(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONTACT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define frictionless dual-mortar surface-to-surface contact. MASTER and SLAVE "
            "must be surface sets. Current slave/master facets are projected onto a common "
            "slave tangent plane and integrated over their physical overlap. Contact uses "
            "augmented-Lagrange normal multipliers and contributes only in NONLINEARSTATIC."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER")   .required()     .doc("Master surface set")
                .key("SLAVE")    .required()     .doc("Slave surface set")
                .key("PENALTY")  .required()     .doc("Normal mortar penalty stiffness")
                .key("CLEARANCE").optional("0")  .doc("Normal clearance subtracted from the signed contact gap")
                .key("FLIP")     .optional("NO") .doc("Flip master surface normals").allowed({"NO", "YES"})
        );

        command.on_enter([&](const fem::io::dsl::Keys& keys) {
            const std::string& master = keys.raw("MASTER");
            const std::string& slave  = keys.raw("SLAVE");
            const fem::Precision k    = keys.get<fem::Precision>("PENALTY");
            const fem::Precision c    = keys.get<fem::Precision>("CLEARANCE");
            const bool flip           = keys.get<bool>("FLIP");

            model.add_contact(master, slave, k, c, flip);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
