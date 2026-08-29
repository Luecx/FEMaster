/**
 * @file register_sfset.cpp
 * @brief Registers part-local and assembly-level surface sets.
 *
 * `SFSET` accepts explicit surface or surface-set references and generated
 * arithmetic ranges. Part/root variants retain sparse surface identifiers in
 * the active Part, while assembly variants are replayed after compilation and
 * resolve optional Instance qualification through compiled surface maps.
 *
 * Scope-dependent behavior is selected declaratively through DSL variants. No
 * separate parser-level assembly flag is required.
 *
 * @see model::Model::add_surfaces_to_part_set
 * @see model::Model::add_surfaces_to_assembly_set
 * @see model::SurfaceRegion
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
 * @brief Registers the `SFSET` command for part-local and compiled assembly scopes.
 */
void register_sfset(dsl::Registry& registry, model::Model& model) {
    registry.command("SFSET", [&](dsl::Command& command) {
        // Allow surface sets in the unqualified root, part and assembly scopes
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            std::string name;
            std::string instance;
        };
        auto ctx = std::make_shared<Context>();

        // Define the set name, optional assembly Instance and generated-range mode
        command.keyword(
            dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                .key("INSTANCE").optional()
                .flag("GENERATE")
        );

        // Store keyword state and create the destination even when no data rows follow
        command.on_enter([&model, ctx](const dsl::ParentInfo& parent, const dsl::Keys& keys) {
            ctx->name     = keys.raw("SFSET");
            ctx->instance = keys.has("INSTANCE") ? keys.raw("INSTANCE") : std::string{};

            if (parent.command == "ASSEMBLY") {
                if (!model._data->compiled) return;

                logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                    "SFSET: instance ", ctx->instance, " is not defined");

                model._data->surface_sets.activate(ctx->name);
                return;
            }

            logging::error(ctx->instance.empty(),
                "SFSET: INSTANCE is only valid at assembly level");
            if (model._data->compiled) return;

            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "SFSET: no active part is available");

            part->surface_sets.activate(ctx->name);
        });

        const auto part_scope     = dsl::Condition::parent_is({"ROOT", "PART"});
        const auto assembly_scope = dsl::Condition::parent_is("ASSEMBLY");
        const auto generated      = dsl::Condition::key_true("GENERATE");
        const auto listed         = dsl::Condition::any_of({
            dsl::Condition::negate(dsl::Condition::key_present("GENERATE")),
            dsl::Condition::key_false("GENERATE")
        });

        // Expand generated ranges directly into sparse part-local surface identifiers
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({part_scope, generated}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("START")
                    .one<ID>().name("END")
                    .one<ID>().name("INC").on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    logging::error(ctx->instance.empty(),
                        "SFSET: INSTANCE is only valid at assembly level");
                    logging::error(inc != 0,
                        "SFSET/GENERATE: increment must not be zero");
                    if (model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SFSET: no active part is available");
                    const auto destination = part->surface_sets.activate(ctx->name);

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
                .range(dsl::LineRange{}.min(0))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("START")
                    .one<ID>().name("END")
                    .one<ID>().name("INC").on_missing(ID{1}).on_empty(ID{1})
                )
                .bind([&model, ctx](ID first, ID last, ID inc) {
                    logging::error(inc != 0,
                        "SFSET/GENERATE: increment must not be zero");
                    if (!model._data->compiled) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "SFSET: instance ", ctx->instance, " is not defined");
                    const auto destination = model._data->surface_sets.activate(ctx->name);

                    for (ID id = first;; id += inc) {
                        const ID value = ctx->instance.empty()
                            ? model.compiled_surface_id(id)
                            : model.compiled_surface_id(id, ctx->instance);
                        destination->add(value);

                        const ID next = static_cast<ID>(id + inc);
                        if (inc > 0 ? next > last : next < last) break;
                    }
                })
            )
        );

        const std::string missing_token = "__MISSING__";

        // Resolve listed part-local surface or surface-set references in the active Part
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
                        "SFSET: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;
                        model.add_surfaces_to_part_set(ctx->name, target);
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
                        "SFSET: instance ", ctx->instance, " is not defined");
                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;

                        const std::string reference = io::reader::qualify_reference(target, ctx->instance);
                        model.add_surfaces_to_assembly_set(ctx->name, reference);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
