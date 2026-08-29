/**
 * @file register_node.cpp
 * @brief Registers part-local and unqualified assembly nodes.
 *
 * Node rows map sparse user identifiers to Cartesian coordinates in the
 * currently active semantic Part. Root-level definitions use the model's
 * default Part, while Instance expansion and dense global enumeration remain
 * deferred to `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers sparse node definitions before model compilation.
 */
void register_node(dsl::Registry& registry, model::Model& model) {
    registry.command("NODE", [&](dsl::Command& command) {
        // Nodes may be defined directly, inside a Part or inside the assembly scope
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        // Every node row is also inserted into the selected node set
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NSET").optional("NALL")
        );

        // Prepare the active node set before processing coordinate rows
        command.on_enter([&model](const dsl::Keys& keys) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "NODE: nodes cannot be added after compile()");

            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "NODE: no active part is available");

            part->node_sets.activate(keys.raw("NSET"));
        });

        // Read sparse node coordinates
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("ID")
                    .one<Precision>().name("X").on_empty(Precision{0}).on_missing(Precision{0})
                    .one<Precision>().name("Y").on_empty(Precision{0}).on_missing(Precision{0})
                    .one<Precision>().name("Z").on_empty(Precision{0}).on_missing(Precision{0})
                )
                .bind([&model](ID id, Precision x, Precision y, Precision z) {
                    model.set_node(id, x, y, z);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
