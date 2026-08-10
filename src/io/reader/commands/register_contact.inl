/**
 * @file register_contact.inl
 * @brief Registers the surface-to-surface dual-mortar `*CONTACT` command.
 *
 * The command resolves one master and one slave surface set and forwards the
 * prescribed normal penalty, clearance and optional master-normal orientation to
 * `Model::add_contact()`. The DSL exposes no geometric search-distance parameter;
 * current master surfaces are tested directly during mortar segmentation.
 *
 * Parsing and validation of generic keyword types remain responsibilities of the
 * DSL engine. Surface-set existence and non-empty regions are validated by the
 * model facade when the command is entered.
 *
 * @see model::Model::add_contact
 * @see constraint::Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

#include <string>

namespace fem::io::reader::commands {

/**
 * Registers the frictionless dual-mortar `*CONTACT` DSL command.
 *
 * The command is valid only at root scope. `MASTER`, `SLAVE` and `PENALTY` are
 * required keyword arguments; `CLEARANCE` defaults to zero and `FLIP` defaults
 * to `NO`. Both regions are interpreted exclusively as surface sets. On command
 * entry the parsed values are forwarded unchanged to `Model::add_contact()`.
 *
 * No data rows are consumed by this command, so it registers one empty variant.
 *
 * @param registry DSL registry receiving the command definition.
 * @param model Model that receives the parsed contact definition.
 */
inline void register_contact(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONTACT", [&](fem::io::dsl::Command& command) {
        // Restrict contact definitions to the model root and describe the
        // surface-to-surface mortar formulation exposed by this command.
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define frictionless dual-mortar surface-to-surface contact. MASTER and SLAVE "
            "must be surface sets. Current slave/master facets are projected onto a common "
            "slave tangent plane and integrated over their physical overlap. Contact uses "
            "augmented-Lagrange normal multipliers and contributes only in NONLINEARSTATIC."
        );

        // Declare the complete keyword interface. PENALTY is the starting AL
        // penalty; later numerical adaptation remains internal to Contact.
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MASTER")   .required()     .doc("Master surface set")
                .key("SLAVE")    .required()     .doc("Slave surface set")
                .key("PENALTY")  .required()     .doc("Normal mortar penalty stiffness")
                .key("CLEARANCE").optional("0")  .doc("Normal clearance subtracted from the signed contact gap")
                .key("FLIP")     .optional("NO") .doc("Flip master surface normals").allowed({"NO", "YES"})
        );

        // Resolve typed keyword values only after the DSL engine has validated
        // required keys, defaults and allowed values.
        command.on_enter([&](const fem::io::dsl::Keys& keys) {
            const std::string&   master = keys.raw("MASTER");
            const std::string&   slave  = keys.raw("SLAVE");
            const fem::Precision k      = keys.get<fem::Precision>("PENALTY");
            const fem::Precision c      = keys.get<fem::Precision>("CLEARANCE");
            const bool           flip   = keys.get<bool>("FLIP");

            model.add_contact(master, slave, k, c, flip);
        });

        // CONTACT has no following data block.
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
