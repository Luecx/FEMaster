/**
 * @file register_nset.inl
 * @brief Registers part-local and assembly-level node sets.
 *
 * `NSET` accepts explicit node or node-set references and generated arithmetic
 * ranges. Part/root variants retain sparse node identifiers in the active Part,
 * while assembly variants are replayed after compilation and resolve optional
 * Instance qualification through compiled node maps.
 *
 * Scope-dependent behavior is selected declaratively through DSL variants. No
 * separate parser-level assembly flag is required.
 *
 * @see model::Model::add_nodes_to_set
 * @see model::NodeRegion
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers the `NSET` command for part-local and compiled assembly scopes.
 */
inline void register_nset(dsl::Registry& registry, model::Model& model) {
    registry.command("NSET", [&](dsl::Command& command) {
        // Allow node sets in the unqualified root, part and assembly scopes
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            std::string name;
            std::string instance;
        };
        auto ctx = std::make_shared<Context>();

        // Define the set name, optional assembly Instance and generated-range mode
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NSET").alternative("NAME").required()
                .key("INSTANCE").optional()
                .flag("GENERATE")
        );

        // Store keyword state shared by all data-layout variants
        command.on_enter([ctx](const dsl::Keys& keys) {
            ctx->name     = keys.raw("NSET");
            ctx->instance = keys.has("INSTANCE") ? keys.raw("INSTANCE") : std::string{};
        });

        const auto part_scope     = dsl::Condition::parent_is({"ROOT", "PART"});
        const auto assembly_scope = dsl::Condition::parent_is("ASSEMBLY");
        const auto generated      = dsl::Condition::key_true("GENERATE");
        const auto listed         = dsl::Condition::any_of({
            dsl::Condition::negate(dsl::Condition::key_present("GENERATE")),
            dsl::Condition::key_false("GENERATE")
        });

        // Expand generated ranges directly into sparse part-local node identifiers
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({part_scope, generated}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("START")
                    .one<ID>().name("END")
                    .one<ID>().name("INC").on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    logging::error(ctx->instance.empty(),
                        "NSET: INSTANCE is only valid at assembly level");
                    logging::error(inc != 0,
                        "NSET/GENERATE: increment must not be zero");
                    if (model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "NSET: no active part is available");
                    const auto destination = part->node_sets.activate(ctx->name);

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
                    .one<ID>().name("START")
                    .one<ID>().name("END")
                    .one<ID>().name("INC").on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    logging::error(inc != 0,
                        "NSET/GENERATE: increment must not be zero");
                    if (!model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "NSET: instance ", ctx->instance, " is not defined");
                    const auto destination = model._data->node_sets.activate(ctx->name);

                    for (ID id = first;; id += inc) {
                        const ID value = ctx->instance.empty()
                            ? model.compiled_node_id(id)
                            : model.compiled_node_id(id, ctx->instance);
                        destination->add(value);

                        const ID next = static_cast<ID>(id + inc);
                        if (inc > 0 ? next > last : next < last) break;
                    }
                })
            )
        );

        const std::string missing_token = "__MISSING__";

        // Resolve listed part-local node or node-set references in the active Part
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::all_of({part_scope, listed}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .fixed<std::string, 32>().name("TARGET")
                        .on_missing(missing_token).on_empty(missing_token)
                )
                .bind([&model, ctx, missing_token](const std::array<std::string, 32>& targets) {
                    logging::error(ctx->instance.empty(),
                        "NSET: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;
                        model.add_nodes_to_set(ctx->name, target);
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
                        .on_missing(missing_token).on_empty(missing_token)
                )
                .bind([&model, ctx, missing_token](const std::array<std::string, 32>& targets) {
                    if (!model._data->compiled) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "NSET: instance ", ctx->instance, " is not defined");

                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;
                        model.add_nodes_to_set(ctx->name, target, ctx->instance);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
