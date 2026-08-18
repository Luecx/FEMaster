/**
 * @file register_surface.inl
 * @brief Registers Abaqus *SURFACE element- and node-based regions.
 *
 * Element-based definitions create FEMaster surface or line entities from an
 * element set or element identifier and an Abaqus side label. Node-based
 * definitions create ordinary FEMaster node regions under the Abaqus surface
 * name. Optional nodal area or weight values are consumed but not stored.
 *
 * @see model::Model::set_surface
 * @see model::NodeRegion
 * @see model::SurfaceRegion
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <charconv>
#include <string>
#include <system_error>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the supported Abaqus `*SURFACE` syntax.
 *
 * `TYPE=ELEMENT` accepts an element set or element identifier followed by a side
 * label such as `S1`, `SPOS` or `SNEG`. `TYPE=NODE` accepts a node set or node
 * identifier and stores the resulting region in the FEMaster node-set registry.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the generated regions.
 */
inline void register_surface(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SURFACE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define an Abaqus element- or node-based surface region.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Abaqus surface name")
                .key("TYPE")
                    .optional("ELEMENT")
                    .allowed({"ELEMENT", "NODE"})
                    .doc("Surface representation")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string& name = keys.raw("NAME");
            const std::string& type = keys.raw("TYPE");

            if (type == "NODE") {
                model._data->node_sets.activate(name);
                return;
            }

            auto surfaces = model._data->surface_sets.activate(name);
            auto lines    = model._data->line_sets.activate(name);
            surfaces->sorted(true).duplicates(false);
            lines->sorted(true).duplicates(false);
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ELEMENT"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Element set or element identifier")
                    .one<std::string>().name("SIDE").desc("Abaqus face label, e.g. S1, SPOS or SNEG")
                )
                .bind([&model](const std::string& target, const std::string& side_token) {
                    int side = 0;
                    if (side_token == "SPOS") {
                        side = 1;
                    } else if (side_token == "SNEG") {
                        side = 2;
                    } else {
                        const std::string numeric =
                            side_token.size() > 1 && side_token[0] == 'S'
                                ? side_token.substr(1)
                                : side_token;
                        const char* begin = numeric.data();
                        const char* end   = begin + numeric.size();
                        const auto [ptr, ec] = std::from_chars(begin, end, side);
                        logging::error(ec == std::errc{} && ptr == end,
                            "SURFACE side '", side_token, "' is not a valid side identifier");
                    }

                    logging::error(side > 0,
                        "SURFACE side must be positive");

                    if (model._data->elem_sets.has(target)) {
                        model.set_surface(target, static_cast<fem::ID>(side));
                        return;
                    }

                    fem::ID element_id{};
                    const char* begin = target.data();
                    const char* end   = begin + target.size();
                    const auto [ptr, ec] = std::from_chars(begin, end, element_id);
                    logging::error(ec == std::errc{} && ptr == end,
                        "SURFACE target '", target, "' is not an element set or element id");
                    model.set_surface(-1, element_id, static_cast<fem::ID>(side));
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"NODE"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or node identifier")
                    .one<fem::Precision>().name("AREA").desc("Ignored Abaqus nodal area/weight")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model](const std::string& target, fem::Precision area) {
                    (void) area;

                    auto destination = model._data->node_sets.get();
                    logging::error(destination != nullptr,
                        "SURFACE TYPE=NODE has no active destination node set");

                    if (model._data->node_sets.has(target)) {
                        auto source = model._data->node_sets.get(target);
                        if (source != destination) {
                            destination->add(*source);
                        }
                        return;
                    }

                    fem::ID node_id{};
                    const char* begin = target.data();
                    const char* end   = begin + target.size();
                    const auto [ptr, ec] = std::from_chars(begin, end, node_id);
                    logging::error(ec == std::errc{} && ptr == end,
                        "SURFACE target '", target, "' is not a node set or node id");
                    destination->add(node_id);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
