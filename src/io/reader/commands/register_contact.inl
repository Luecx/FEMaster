/**
 * @file register_contact.inl
 * @brief Registers frictionless node-to-surface penalty contact.
 *
 * `MASTER` and `SLAVE` resolve to surface sets. The nodes contained in the slave
 * surface set are used as contact nodes, while every face in the master surface
 * set is tested directly as a potential master face. No search-distance,
 * surface-to-surface, mortar or augmented-Lagrange formulation is exposed.
 *
 * @see model::Model::add_contact
 * @see constraint::Contact
 *
 * @author Finn Eggers
 * @date 11.08.2026
 */

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

#include <string>

namespace fem::io::reader::commands {

/**
 * Registers the frictionless node-to-surface `*CONTACT` command.
 *
 * `MASTER`, `SLAVE` and `PENALTY` are required. `CLEARANCE` defaults to zero and
 * `FLIP` optionally reverses the master face normal. Contact contributes only in
 * nonlinear static assembly.
 *
 * @param registry DSL registry receiving the command definition.
 * @param model Model that receives the parsed contact definition.
 */
inline void register_contact(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONTACT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define frictionless node-to-surface penalty contact. MASTER and SLAVE must be "
            "surface sets. Every unique slave-surface node is projected onto all master "
            "faces and the closest bounded master-face projection is used. Contact acts "
            "only for negative normal gap and contributes only in NONLINEARSTATIC."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER")   .required()     .doc("Master surface set")
                .key("SLAVE")    .required()     .doc("Surface set providing slave contact nodes")
                .key("PENALTY")  .required()     .doc("Normal penalty stiffness per unit area")
                .key("CLEARANCE").optional("0")  .doc("Clearance subtracted from the signed normal gap")
                .key("FLIP")     .optional("NO") .doc("Flip master surface normals").allowed({"NO", "YES"})
        );

        command.on_enter([&](const fem::io::dsl::Keys& keys) {
            const std::string&   master = keys.raw("MASTER");
            const std::string&   slave  = keys.raw("SLAVE");
            const fem::Precision k      = keys.get<fem::Precision>("PENALTY");
            const fem::Precision c      = keys.get<fem::Precision>("CLEARANCE");
            const bool           flip   = keys.get<bool>("FLIP");

            model.add_contact(master, slave, k, c, flip);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
