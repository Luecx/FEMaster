/**
 * @file register_surface.cpp
 * @brief Registers part-local and assembly-level node- and element-based surfaces.
 *
 * The command provides the shared `SURFACE` grammar used by FEMaster and Abaqus
 * input decks. Element-based rows associate an element or element set with one
 * of its boundary sides and construct the corresponding surface or line region.
 * Node-based rows populate a named `NodeRegion`.
 *
 * Part/root and assembly semantics are represented by separate DSL variants.
 * Part-local definitions retain sparse semantic identifiers before compilation,
 * while assembly variants are replayed afterwards and resolve optional Instance
 * qualification through the compiled local-to-global maps.
 *
 * @see model::NodeRegion
 * @see model::SurfaceRegion
 * @see model::LineRegion
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

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

namespace dsl = fem::io::dsl;

/**
 * @brief Registers node- and element-based `SURFACE` definitions.
 */
void register_surface(dsl::Registry& registry, model::Model& model) {
    registry.command("SURFACE", [&](dsl::Command& command) {
        // Allow surfaces in the unqualified root, part and assembly scopes
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));

        struct Context {
            std::string name;
            std::string type;
            std::string instance;
        };
        auto ctx = std::make_shared<Context>();

        // Define naming, surface representation and optional assembly qualification
        command.keyword(
            dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                .key("TYPE").optional("ELEMENT").allowed({"ELEMENT", "NODE"})
                .key("INSTANCE").optional()
        );

        // Store keyword state shared by all surface-layout variants
        command.on_enter([ctx](const dsl::Keys& keys) {
            ctx->name     = keys.raw("SFSET");
            ctx->type     = keys.raw("TYPE");
            ctx->instance = keys.has("INSTANCE") ? keys.raw("INSTANCE") : std::string{};
        });

        const auto part_scope      = dsl::Condition::parent_is({"ROOT", "PART"});
        const auto assembly_scope  = dsl::Condition::parent_is("ASSEMBLY");
        const auto element_surface = dsl::Condition::key_equals("TYPE", {"ELEMENT"});
        const auto node_surface    = dsl::Condition::key_equals("TYPE", {"NODE"});

        // Translate Abaqus/FEMaster side labels into the element boundary index
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

        // Materialize one compiled element boundary as a surface or line entity
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

        // Preserve the FEMaster explicit surface-id form for part-local definitions
        command.variant(dsl::Variant::make()
            .rank(20)
            .when(dsl::Condition::all_of({part_scope, element_surface}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<ID>().name("ID")
                    .one<ID>().name("ELEM_ID")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side](ID id, ID element_id, const std::string& side_token) {
                    logging::error(ctx->instance.empty(),
                        "SURFACE: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SURFACE: no active part is available");
                    part->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                    part->line_sets.activate(ctx->name)->sorted(true).duplicates(false);

                    model.set_surface(id, element_id, parse_side(side_token));
                })
            )
        );

        // Resolve part-local element or element-set references against sparse topology
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({part_scope, element_surface}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side](const std::string& target, const std::string& side_token) {
                    logging::error(ctx->instance.empty(),
                        "SURFACE: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    const ID side = static_cast<ID>(parse_side(side_token));
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SURFACE: no active part is available");

                    part->surface_sets.activate(ctx->name)->sorted(true).duplicates(false);
                    part->line_sets.activate(ctx->name)->sorted(true).duplicates(false);

                    if (part->elem_sets.has(target)) {
                        model.set_surface(target, side);
                    } else {
                        model.set_surface(-1, io::reader::parse_local_id(target, "SURFACE"), side);
                    }
                })
            )
        );

        // Resolve assembly element references and materialize their compiled boundaries
        command.variant(dsl::Variant::make()
            .rank(10)
            .when(dsl::Condition::all_of({assembly_scope, element_surface}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<std::string>().name("SIDE")
                )
                .bind([&model, ctx, parse_side, add_compiled_boundary](const std::string& target,
                                                                       const std::string& side_token) {
                    if (!model._data->compiled) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "SURFACE: instance ", ctx->instance, " is not defined");

                    const ID side = static_cast<ID>(parse_side(side_token));
                    const std::string reference = io::reader::qualify_reference(target, ctx->instance);

                    if (model._data->elem_sets.has(reference)) {
                        const auto elements = model._data->elem_sets.get(reference);
                        logging::error(elements != nullptr,
                            "SURFACE: element set ", reference, " is not initialized");

                        for (const ID element_id : *elements) {
                            add_compiled_boundary(element_id, side, ctx->name);
                        }
                    } else {
                        add_compiled_boundary(model.compiled_element_id(reference), side, ctx->name);
                    }
                })
            )
        );

        // Populate a part-local node-based surface through sparse node references
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::all_of({part_scope, node_surface}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<Precision>().name("AREA").on_missing(Precision{1}).on_empty(Precision{1})
                )
                .bind([&model, ctx](const std::string& target, Precision area) {
                    (void) area;
                    logging::error(ctx->instance.empty(),
                        "SURFACE: INSTANCE is only valid at assembly level");
                    if (model._data->compiled) return;

                    model.add_nodes_to_part_set(ctx->name, target);
                })
            )
        );

        // Populate an assembly node-based surface through compiled node references
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::all_of({assembly_scope, node_surface}))
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<Precision>().name("AREA").on_missing(Precision{1}).on_empty(Precision{1})
                )
                .bind([&model, ctx](const std::string& target, Precision area) {
                    (void) area;
                    if (!model._data->compiled) return;

                    logging::error(ctx->instance.empty() || model._data->instances.has(ctx->instance),
                        "SURFACE: instance ", ctx->instance, " is not defined");

                    const std::string reference = io::reader::qualify_reference(target, ctx->instance);
                    model.add_nodes_to_assembly_set(ctx->name, reference);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
