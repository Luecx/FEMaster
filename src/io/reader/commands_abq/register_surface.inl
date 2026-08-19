/**
 * @file register_surface.inl
 * @brief Registers Abaqus part-local and assembly-level surfaces.
 *
 * Abaqus element-based surface entries are translated from element-set and side
 * labels into FEMaster surface identifiers and `SurfaceRegion` objects. Before
 * compilation, definitions remain part-local; assembly definitions are replayed
 * later so Instance-qualified references can resolve through compiled maps.
 *
 * The command preserves named surface groups used by pressure, traction,
 * contact, tie and coupling translations without duplicating element geometry.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <charconv>
#include <memory>
#include <string>
#include <system_error>

#include "../reference.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_surface(fem::io::dsl::Registry& registry,
                             model::Model& model,
                             std::shared_ptr<bool> assembly_scope) {
    registry.command("SURFACE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            bool assembly = false;
            std::string name;
            std::string type;
            model::NodeRegion::Ptr node_destination = nullptr;
        };
        auto ctx = std::make_shared<Context>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("TYPE").optional("ELEMENT").allowed({"ELEMENT", "NODE"})
        );
        command.on_enter([&model, assembly_scope, ctx](const fem::io::dsl::Keys& keys) {
            ctx->assembly         = *assembly_scope;
            ctx->name             = keys.raw("NAME");
            ctx->type             = keys.raw("TYPE");
            ctx->node_destination = nullptr;

            if (ctx->assembly) {
                if (!model._data->compiled) return;
                if (ctx->type == "NODE") {
                    ctx->node_destination = model._data->node_sets.activate(ctx->name);
                } else {
                    model._data->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                    model._data->line_sets.activate(ctx->name)->sorted(true).duplicates(false);
                }
                return;
            }

            if (model._data->compiled) return;
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "SURFACE: no active part is available");
            if (ctx->type == "NODE") {
                ctx->node_destination = part->node_sets.activate(ctx->name);
            } else {
                part->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                part->line_sets.activate(ctx->name)->sorted(true).duplicates(false);
            }
        });

        auto parse_side = [](const std::string& side_token) -> ID {
            if (side_token == "SPOS") return 1;
            if (side_token == "SNEG") return 2;
            const std::string numeric = side_token.size() > 1 && side_token[0] == 'S'
                ? side_token.substr(1) : side_token;
            int side{};
            const char* begin = numeric.data();
            const char* end   = begin + numeric.size();
            const auto [ptr, ec] = std::from_chars(begin, end, side);
            logging::error(ec == std::errc{} && ptr == end && side > 0,
                "SURFACE: side ", side_token, " is not a valid positive side identifier");
            return static_cast<ID>(side);
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
                "SURFACE: side ", side, " of compiled element ", element_id,
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
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ELEMENT"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side, add_compiled_boundary](const std::string& target,
                                                                       const std::string& side_token) {
                    const ID side = parse_side(side_token);
                    if (ctx->assembly) {
                        if (!model._data->compiled) return;
                        if (model._data->elem_sets.has(target)) {
                            const auto elements = model._data->elem_sets.get(target);
                            logging::error(elements != nullptr,
                                "SURFACE: element set ", target, " is not initialized");
                            for (const ID element_id : *elements) add_compiled_boundary(element_id, side, ctx->name);
                        } else {
                            add_compiled_boundary(io::reader::compiled_element_id(model, target), side, ctx->name);
                        }
                        return;
                    }

                    if (model._data->compiled) return;
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SURFACE: no active part is available");
                    if (part->elem_sets.has(target)) model.set_surface(target, side);
                    else model.set_surface(-1, io::reader::parse_local_id(target, "SURFACE"), side);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"NODE"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<Precision>().name("AREA").on_missing(Precision{1}).on_empty(Precision{1})
                )
                .bind([&model, ctx](const std::string& target, Precision area) {
                    (void) area;
                    if (!ctx->node_destination) return;
                    if (ctx->assembly) {
                        io::reader::add_compiled_reference(
                            model._data->node_sets,
                            ctx->node_destination,
                            target,
                            "",
                            [&model](const std::string& reference) {
                                return io::reader::compiled_node_id(model, reference);
                            }
                        );
                        return;
                    }

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SURFACE: no active part is available");
                    if (part->node_sets.has(target)) {
                        const auto source = part->node_sets.get(target);
                        logging::error(source != nullptr,
                            "SURFACE: node set ", target, " is not initialized");
                        if (source != ctx->node_destination) ctx->node_destination->add(*source);
                    } else {
                        ctx->node_destination->add(io::reader::parse_local_id(target, "SURFACE"));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
