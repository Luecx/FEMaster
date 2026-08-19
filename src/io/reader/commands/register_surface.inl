/**
 * @file register_surface.inl
 * @brief Registers part-local and post-compile assembly surfaces.
 *
 * Surface rows associate a sparse surface identifier with an element and one
 * of its boundary sides. The command accepts numeric and conventional side
 * labels, creates the element-specific surface representation and stores it in
 * the active Part before compilation.
 *
 * Assembly-level definitions are replayed after compilation so Instance-qualified
 * element references can be translated to dense global IDs without duplicating
 * semantic topology.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <charconv>
#include <cctype>
#include <memory>
#include <string>
#include <system_error>

#include "../reference.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_surface(fem::io::dsl::Registry& registry,
                             model::Model& model,
                             std::shared_ptr<bool> assembly_scope) {
    registry.command("SURFACE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            bool assembly = false;
            std::string name;
            std::string instance;
        };
        auto ctx = std::make_shared<Context>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                .key("INSTANCE").optional()
        );
        command.on_enter([&model, assembly_scope, ctx](const fem::io::dsl::Keys& keys) {
            ctx->assembly = *assembly_scope;
            ctx->name     = keys.raw("SFSET");
            ctx->instance = keys.has("INSTANCE") ? keys.raw("INSTANCE") : std::string{};

            logging::error(ctx->assembly || ctx->instance.empty(),
                "SURFACE: INSTANCE is only valid at assembly level");

            if (ctx->assembly) {
                if (!model._data->compiled) return;
                logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                    "SURFACE: instance ", ctx->instance, " is not defined");
                model._data->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                model._data->line_sets.activate(ctx->name)->sorted(true).duplicates(false);
            } else if (!model._data->compiled) {
                const auto part = model._data->parts.get();
                logging::error(part != nullptr,
                    "SURFACE: no active part is available");
                part->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                part->line_sets.activate(ctx->name)->sorted(true).duplicates(false);
            }
        });

        auto parse_side = [](const std::string& side_token) -> int {
            if (side_token.empty()) return 1;

            std::string token = side_token;
            for (char& c : token) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
            if (token == "SPOS") return 1;
            if (token == "SNEG") return 2;

            const std::string numeric = token[0] == 'S' && token.size() > 1 ? token.substr(1) : token;
            int side{};
            const char* begin = numeric.data();
            const char* end   = begin + numeric.size();
            const auto [ptr, ec] = std::from_chars(begin, end, side);
            logging::error(ec == std::errc{} && ptr == end,
                "SURFACE: side ", side_token, " is not a valid face identifier");
            return side;
        };

        auto add_compiled_boundary = [&model](ID element_id, ID side, const std::string& name) {
            logging::error(element_id >= 0 && static_cast<std::size_t>(element_id) < model._data->elements.size(),
                "SURFACE: compiled element ", element_id, " is out of range");
            logging::error(model._data->elements[static_cast<std::size_t>(element_id)] != nullptr,
                "SURFACE: compiled element ", element_id, " is not configured");

            const auto element = model._data->elements[static_cast<std::size_t>(element_id)];
            auto surface = element->surface(side);
            auto line    = element->line(side);
            logging::error(surface != nullptr || line != nullptr,
                "SURFACE: boundary ", side, " of compiled element ", element_id,
                " provides neither a surface nor a line");

            if (surface) {
                const ID surface_id = static_cast<ID>(model._data->surfaces.size());
                model._data->surfaces.push_back(std::move(surface));
                model._data->surface_sets.activate(name)->add(surface_id);
            }
            if (line) {
                const ID line_id = static_cast<ID>(model._data->lines.size());
                model._data->lines.push_back(std::move(line));
                model._data->line_sets.activate(name)->add(line_id);
            }
        };

        command.variant(fem::io::dsl::Variant::make()
            .rank(20)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<ID>().name("ID")
                    .one<ID>().name("ELEM_ID")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side](ID id, ID element_id, const std::string& side_token) {
                    logging::error(!ctx->assembly,
                        "SURFACE: explicit surface ids are not supported at assembly level");
                    if (!model._data->compiled) {
                        model.set_surface(id, element_id, parse_side(side_token));
                    }
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .rank(10)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side, add_compiled_boundary](const std::string& target,
                                                                       const std::string& side_token) {
                    const ID side = static_cast<ID>(parse_side(side_token));
                    if (ctx->assembly) {
                        if (!model._data->compiled) return;

                        const std::string reference = io::reader::qualify_reference(target, ctx->instance);
                        if (model._data->elem_sets.has(reference)) {
                            const auto elements = model._data->elem_sets.get(reference);
                            logging::error(elements != nullptr,
                                "SURFACE: element set ", reference, " is not initialized");
                            for (const ID element_id : *elements) {
                                add_compiled_boundary(element_id, side, ctx->name);
                            }
                        } else {
                            add_compiled_boundary(io::reader::compiled_element_id(model, reference), side, ctx->name);
                        }
                        return;
                    }

                    if (model._data->compiled) return;
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SURFACE: no active part is available");
                    if (part->elem_sets.has(target)) {
                        model.set_surface(target, side);
                    } else {
                        model.set_surface(-1, io::reader::parse_local_id(target, "SURFACE"), side);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
