/**
 * @file register_elset.cpp
 * @brief Registers the ELSET input command in semantic Part and compiled Assembly scopes.
 *
 * `ELSET` creates named `model::ElementRegion` objects from either explicit
 * references or inclusive arithmetic ranges. An explicit token may identify one
 * element or another element set in the same identifier space. `GENERATE` rows
 * provide a start identifier, end identifier and optional non-zero increment.
 *
 * Root- and Part-level definitions are processed before model compilation and
 * retain sparse Part-local element identifiers. Assembly definitions are
 * replayed after compilation; scalar references are mapped through the selected
 * Instance, while named source sets already contain dense assembly identifiers.
 * The Reader applies a separate `INSTANCE` keyword to otherwise unqualified
 * tokens, and `model::Model` performs set composition and identifier mapping.
 *
 * @see model::Model::add_elements_to_part_set
 * @see model::Model::add_elements_to_assembly_set
 * @see model::ElementRegion
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers explicit and generated element-set definitions.
 *
 * Command entry stores the destination name and optional Instance and creates
 * the destination region in the identifier space selected by the parent scope.
 * Four DSL variants distinguish listed and generated rows for semantic Part and
 * compiled Assembly processing. Creating the destination during command entry
 * preserves an explicitly declared empty set.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model providing the active Part, compiled Instance maps and
 *              element-set registries.
 */
void register_elset(dsl::Registry& registry, model::Model& model) {
    registry.command("ELSET", [&](dsl::Command& command) {
        // Element sets may address semantic topology or the compiled assembly
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        // Expose the two supported row forms and their scope-dependent semantics
        command.doc(
            "Create an element set from explicit element or element-set references, "
            "or from inclusive START, END, INC ranges when GENERATE is present."
        );

        /**
         * @brief Keyword state shared by the scope- and row-specific callbacks.
         *
         * `name` identifies the destination region. `instance` is empty for
         * Root/Part definitions and optionally qualifies Assembly references.
         */
        struct Context {
            // Destination element-set name
            std::string name;

            // Optional Instance qualifier for Assembly rows
            std::string instance;
        };
        auto ctx = std::make_shared<Context>();

        // Define the set name, optional assembly Instance and generated-range mode
        command.keyword(
            dsl::KeywordSpec::make()
                .key("ELSET").alternative("NAME").required().doc("Destination element-set name")
                .key("INSTANCE").optional().doc("Instance qualifying unqualified Assembly references")
                .flag("GENERATE").doc("Interpret each row as an inclusive arithmetic range")
        );

        // Store command-wide syntax state and materialize the destination even
        // when the command contains no listed data rows.
        command.on_enter([&model, ctx](const dsl::ParentInfo& parent, const dsl::Keys& keys) {
            ctx->name     = keys.raw("ELSET");
            ctx->instance = keys.raw("INSTANCE");

            if (parent.command == "ASSEMBLY") {
                // Assembly commands are replayed after topology compilation
                if (!model._data->compiled) return;

                logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                    "ELSET: instance ", ctx->instance, " is not defined");

                // Preserve an explicitly declared empty Assembly element set
                model._data->elem_sets.activate(ctx->name);
                return;
            }

            // Root and Part definitions cannot carry Assembly qualification
            logging::error(ctx->instance.empty(),
                "ELSET: INSTANCE is only valid at assembly level");
            if (model._data->compiled) return;

            // Preserve an explicitly declared empty set in the active Part
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "ELSET: no active part is available");

            part->elem_sets.activate(ctx->name);
        });

        // Select one of four variants from the parent scope and GENERATE flag
        const auto part_scope     = dsl::Condition::parent_is({"ROOT", "PART"});
        const auto assembly_scope = dsl::Condition::parent_is("ASSEMBLY");
        const auto generated      = dsl::Condition::key_true("GENERATE");
        const auto listed         = dsl::Condition::any_of({
            dsl::Condition::negate(dsl::Condition::key_present("GENERATE")),
            dsl::Condition::key_false("GENERATE")
        });

        // Expand generated ranges directly into sparse part-local element identifiers
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({part_scope, generated}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("START").desc("First Part-local element identifier")
                    .one<ID>().name("END"  ).desc("Last Part-local element identifier")
                    .one<ID>().name("INC"  ).desc("Non-zero identifier increment")
                        .on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    // Validate the Part scope and arithmetic progression
                    logging::error(ctx->instance.empty(),
                        "ELSET: INSTANCE is only valid at assembly level");
                    logging::error(inc != 0,
                        "ELSET/GENERATE: increment must not be zero");
                    if (model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "ELSET: no active part is available");
                    const auto destination = part->elem_sets.activate(ctx->name);

                    // Insert the inclusive progression as sparse Part-local identifiers
                    for (ID id = first;; id += inc) {
                        destination->add(id);

                        const ID next = static_cast<ID>(id + inc);
                        if (inc > 0 ? next > last : next < last) break;
                    }
                })
            )
        );

        // Expand generated assembly ranges through the compiled Instance mapping
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({assembly_scope, generated}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("START").desc("First Instance-local element identifier")
                    .one<ID>().name("END"  ).desc("Last Instance-local element identifier")
                    .one<ID>().name("INC"  ).desc("Non-zero identifier increment")
                        .on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    // Validate the compiled scope and arithmetic progression
                    logging::error(inc != 0,
                        "ELSET/GENERATE: increment must not be zero");
                    if (!model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "ELSET: instance ", ctx->instance, " is not defined");
                    const auto destination = model._data->elem_sets.activate(ctx->name);

                    // Map every Instance-local identifier into dense assembly space
                    for (ID id = first;; id += inc) {
                        const ID compiled_id = ctx->instance.empty()
                            ? model.compiled_element_id(id)
                            : model.compiled_element_id(id, ctx->instance);
                        destination->add(compiled_id);

                        const ID next = static_cast<ID>(id + inc);
                        if (inc > 0 ? next > last : next < last) break;
                    }
                })
            )
        );

        // Fixed-width listed rows use a sentinel for absent trailing columns
        const std::string missing_token = "__MISSING__";

        // Resolve listed part-local element or element-set references in the active Part
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::all_of({part_scope, listed}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .fixed<std::string, 32>().name("TARGET")
                        .desc("Part-local element identifier or element-set name")
                        .on_missing(missing_token).on_empty(missing_token)
                )
                .bind([&model, ctx, missing_token](const std::array<std::string, 32>& targets) {
                    // Process this variant only in semantic Part space
                    logging::error(ctx->instance.empty(),
                        "ELSET: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    // Delegate named-set composition and scalar parsing to Model
                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;
                        model.add_elements_to_part_set(ctx->name, target);
                    }
                })
            )
        );

        // Resolve listed assembly references with optional Instance qualification
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::all_of({assembly_scope, listed}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .fixed<std::string, 32>().name("TARGET")
                        .desc("Assembly element reference or element-set name")
                        .on_missing(missing_token).on_empty(missing_token)
                )
                .bind([&model, ctx, missing_token](const std::array<std::string, 32>& targets) {
                    // Process this variant only after assembly compilation
                    if (!model._data->compiled) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "ELSET: instance ", ctx->instance, " is not defined");

                    // Apply the keyword-level Instance only to unqualified source tokens
                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;

                        const std::string reference = io::reader::qualify_reference(target, ctx->instance);
                        model.add_elements_to_assembly_set(ctx->name, reference);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
