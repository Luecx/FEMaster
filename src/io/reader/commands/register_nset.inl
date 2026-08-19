/**
 * @file register_nset.inl
 * @brief Registers part-local and assembly-level node sets.
 *
 * `NSET` accepts explicit identifiers and generated arithmetic ranges. Before
 * compilation it stores sparse local node IDs in the active Part; during the
 * assembly pass it resolves optional Instance-qualified references through the
 * compiled local-to-global maps.
 *
 * The resulting `NodeRegion` is registered under its deck name for supports,
 * loads, couplings and later assembly commands.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_nset(fem::io::dsl::Registry& registry,
                          model::Model& model,
                          std::shared_ptr<bool> assembly_scope) {
    registry.command("NSET", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            bool assembly = false;
            std::string instance;
            model::NodeRegion::Ptr destination = nullptr;
        };
        auto ctx = std::make_shared<Context>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").alternative("NAME").required()
                .key("INSTANCE").optional()
                .flag("GENERATE")
        );
        command.on_enter([&model, assembly_scope, ctx](const fem::io::dsl::Keys& keys) {
            ctx->assembly    = *assembly_scope;
            ctx->instance    = keys.has("INSTANCE") ? keys.raw("INSTANCE") : std::string{};
            ctx->destination = nullptr;

            logging::error(ctx->assembly || ctx->instance.empty(),
                "NSET: INSTANCE is only valid at assembly level");

            if (ctx->assembly) {
                if (!model._data->compiled) return;
                logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                    "NSET: instance ", ctx->instance, " is not defined");
                ctx->destination = model._data->node_sets.activate(keys.raw("NSET"));
            } else if (!model._data->compiled) {
                const auto part = model._data->parts.get();
                logging::error(part != nullptr,
                    "NSET: no active part is available");
                ctx->destination = part->node_sets.activate(keys.raw("NSET"));
            }
        });

        const std::string missing_token = "__MISSING__";

        command.variant(fem::io::dsl::Variant::make()
            .rank(10)
            .when(fem::io::dsl::Condition::key_true("GENERATE"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::ID>().name("START")
                    .one<fem::ID>().name("END")
                    .one<fem::ID>().name("INC").on_missing(fem::ID{1}).on_empty(fem::ID{1})
                )
                .bind([&model, ctx](fem::ID first, fem::ID last, fem::ID inc) {
                    logging::error(inc != 0,
                        "NSET/GENERATE: increment must not be zero");
                    if (!ctx->destination) return;
                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) return;

                    for (fem::ID id = first;; id += inc) {
                        const fem::ID value = ctx->assembly
                            ? (ctx->instance.empty() ? model.compiled_node_id(id)
                                                     : model.compiled_node_id(id, ctx->instance))
                            : id;
                        ctx->destination->add(value);

                        const fem::ID next = static_cast<fem::ID>(id + inc);
                        if (inc > 0 ? next > last : next < last) break;
                    }
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::any_of({
                fem::io::dsl::Condition::negate(fem::io::dsl::Condition::key_present("GENERATE")),
                fem::io::dsl::Condition::key_false("GENERATE")
            }))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 32>().name("TARGET")
                        .on_missing(missing_token).on_empty(missing_token)
                )
                .bind([&model, ctx, missing_token](const std::array<std::string, 32>& targets) {
                    if (!ctx->destination) return;
                    for (const std::string& target : targets) {
                        if (target == missing_token) continue;
                        if (ctx->assembly) {
                            io::reader::add_compiled_reference(
                                model._data->node_sets,
                                ctx->destination,
                                target,
                                ctx->instance,
                                [&model](const std::string& reference) {
                                    return io::reader::compiled_node_id(model, reference);
                                }
                            );
                        } else {
                            ctx->destination->add(io::reader::parse_local_id(target, "NSET"));
                        }
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
