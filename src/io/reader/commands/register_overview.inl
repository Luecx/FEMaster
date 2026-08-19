/**
 * @file register_overview.inl
 * @brief Registers the command that prints the compiled model overview.
 *
 * The root-level `OVERVIEW` command delegates diagnostic output to `Model`
 * after the parser has compiled part and instance topology into assembly data.
 *
 * The report can therefore compare retained semantic topology with dense model
 * fields, regions and solver-facing capacities. The command is read-only and
 * does not trigger compilation or mutate model state.
 *
 * @see fem::model::Model
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

/**
 * Registers the root-level command that emits the model diagnostic overview.
 *
 * The parser executes this command during the post-compile data pass, so the
 * member function can report both retained semantic topology and the resulting
 * dense assembly representation.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model inspected when an `OVERVIEW` command is encountered.
 */
inline void register_overview(fem::io::dsl::Registry& registry, fem::model::Model& model) {
    registry.command("OVERVIEW", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Print the semantic topology, compiled assembly, definitions, "
            "constraints and boundary-condition collectors of the current model."
        );

        command.on_enter([&model](const fem::io::dsl::Keys&) {
            model.print_overview();
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
