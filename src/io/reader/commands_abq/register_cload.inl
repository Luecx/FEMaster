/**
 * @file register_cload.inl
 * @brief Registers Abaqus *CLOAD concentrated nodal load definitions.
 *
 * Concentrated loads are translated directly into the FEMaster load collector of
 * the single supported Abaqus analysis step. Node-set targets are expanded while
 * parsing so nodal `*TRANSFORM` coordinate systems can be assigned to each
 * generated load.
 *
 * @see ParserAbqState
 * @see model::Model::add_cload
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <charconv>
#include <memory>
#include <string>
#include <system_error>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the supported Abaqus concentrated-load syntax.
 *
 * Each data line uses `node-or-nset, dof, magnitude`. Structural DOFs 1 through 6
 * map to the three force and three moment components. Named amplitudes and the
 * step-level default amplitude are resolved according to the active procedure.
 * Follower concentrated loads and imaginary harmonic components are unsupported.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active step and nodal transforms.
 */
inline void register_cload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("CLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define concentrated nodal loads in the active Abaqus step.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .flag("FOLLOWER").doc("Unsupported follower concentrated-load flag")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && parser.active_loadcase(),
                "CLOAD must appear after a supported procedure inside STEP");
            logging::error(parser.active_loadcase_type() != "EIGENFREQ",
                "CLOAD is not supported in a FREQUENCY step");
            logging::error(!keys.has("FOLLOWER"),
                "CLOAD FOLLOWER is not supported");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "CLOAD REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "CLOAD IMAGINARY is not supported by the real-load harmonic solver");

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or node identifier")
                    .one<int>().name("DOF").desc("Abaqus generalized degree of freedom 1--6")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Load magnitude")
                )
                .bind([&parser, amplitude](const std::string& target,
                                           int dof,
                                           fem::Precision magnitude) {
                    logging::error(dof >= 1 && dof <= 6,
                        "CLOAD supports only structural DOFs 1 through 6");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    magnitude *= scale;
                    if (magnitude == fem::Precision(0)) {
                        return;
                    }

                    fem::Vec6 load = fem::Vec6::Zero();
                    load[dof - 1] = magnitude;

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    const auto add_to_node = [&](fem::ID node_id) {
                        std::string orientation;
                        const auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            orientation = transform->second;
                        }
                        model.add_cload(node_id, load, orientation, resolved_amplitude);
                    };

                    if (model._data->node_sets.has(target)) {
                        for (const fem::ID node_id : *model._data->node_sets.get(target)) {
                            add_to_node(node_id);
                        }
                        return;
                    }

                    fem::ID node_id{};
                    const char* begin = target.data();
                    const char* end   = begin + target.size();
                    const auto [ptr, ec] = std::from_chars(begin, end, node_id);
                    logging::error(ec == std::errc{} && ptr == end,
                        "CLOAD target '", target, "' is not a node set or node id");
                    add_to_node(node_id);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
