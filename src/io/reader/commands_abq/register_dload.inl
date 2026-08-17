/**
 * @file register_dload.inl
 * @brief Registers supported step-dependent Abaqus *DLOAD definitions.
 *
 * The current element-based distributed-load subset supports `GRAV`. Definitions
 * are retained logically across steps so `OP=MOD` can update an existing load
 * and `OP=NEW` can replace the complete active DLOAD category without carrying
 * any mechanical solver state between FEMaster load cases.
 *
 * Abaqus identifies a gravity load by target, load type and specified load
 * direction. Repeated definitions of that same identity inside one step are
 * stored separately and add during materialization; an ambiguous propagated
 * multi-definition must be replaced with `OP=NEW` rather than modified.
 *
 * @see ParserAbqState
 * @see ParserAbqDLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus gravity-load history records.
 *
 * The logical identity consists of the original element/element-set target,
 * `GRAV` type and the three supplied gravity-vector components. One propagated
 * definition can be redefined with `OP=MOD`; repeated definitions in the current
 * step remain separate additive load conditions. An empty `OP=NEW` block removes
 * every active DLOAD.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining active distributed-load definitions.
 */
inline void register_dload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define or modify active Abaqus GRAV distributed loads.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("MOD preserves active DLOADs; NEW replaces the complete DLOAD category")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("DLOAD must appear after a supported procedure inside STEP");
            }
            if (state.procedure == "EIGENFREQ") {
                throw std::runtime_error("DLOAD is not supported in a FREQUENCY step");
            }
            if (keys.has("REAL") && keys.has("IMAGINARY")) {
                throw std::runtime_error("DLOAD REAL and IMAGINARY are mutually exclusive");
            }
            if (keys.has("IMAGINARY")) {
                throw std::runtime_error("DLOAD IMAGINARY is not supported by the real-load harmonic solver");
            }

            const std::string op = keys.raw("OP");
            if (!state.dload_op.empty() && state.dload_op != op) {
                throw std::runtime_error("All DLOAD cards in one STEP must use the same OP value");
            }
            if (state.dload_op.empty()) {
                state.dload_op = op;
                if (op == "NEW") {
                    state.dloads.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            if (!amplitude->empty() && !parser.model()._data->amplitudes.has(*amplitude)) {
                throw std::runtime_error("DLOAD references unknown amplitude '" + *amplitude + "'");
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Element set, element id or blank for all elements")
                        .on_empty(std::string{"EALL"})
                    .one<std::string>().name("TYPE").desc("Supported load label GRAV")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Gravity magnitude")
                    .fixed<fem::Precision, 3>().name("DIRECTION").desc("Gravity vector components")
                )
                .bind([&parser, amplitude](const std::string& target,
                                           const std::string& type,
                                           fem::Precision magnitude,
                                           const std::array<fem::Precision, 3>& direction) {
                    if (type != "GRAV") {
                        throw std::runtime_error(
                            "DLOAD currently supports only GRAV; use SURFACE + DSLOAD for surface pressure/traction"
                        );
                    }

                    auto& state = parser.abaqus_state();
                    auto first = state.dloads.end();
                    int matches = 0;
                    bool defined_this_step = false;

                    for (auto it = state.dloads.begin(); it != state.dloads.end(); ++it) {
                        if (it->target != target || it->type != type || it->direction != direction) {
                            continue;
                        }
                        if (first == state.dloads.end()) {
                            first = it;
                        }
                        ++matches;
                        defined_this_step = defined_this_step || it->modified_step == state.step_index;
                    }

                    if (defined_this_step || matches == 0) {
                        state.dloads.push_back(ParserAbqDLoad{
                            target,
                            type,
                            magnitude,
                            Precision(0),
                            direction,
                            *amplitude,
                            state.step_index
                        });
                        return;
                    }
                    if (matches > 1) {
                        throw std::runtime_error(
                            "Multiple active DLOADs use the same target, type and direction; use OP=NEW to redefine them"
                        );
                    }

                    first->previous_magnitude = first->magnitude;
                    first->magnitude          = magnitude;
                    first->amplitude          = *amplitude;
                    first->modified_step      = state.step_index;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
