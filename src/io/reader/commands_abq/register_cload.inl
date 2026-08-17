/**
 * @file register_cload.inl
 * @brief Registers step-dependent Abaqus *CLOAD definitions.
 *
 * Concentrated loads are retained as logical Abaqus records across steps rather
 * than being converted immediately to FEMaster load objects. `OP=MOD` updates or
 * adds target/DOF records; `OP=NEW` clears all active concentrated loads before
 * the current step supplies its replacement definitions.
 *
 * Deferring materialization preserves the original node/node-set identity needed
 * by Abaqus modification semantics. At `*END STEP` node sets are expanded and
 * any nodal `*TRANSFORM` is copied to the concrete FEMaster `CLoad` objects.
 *
 * @see ParserAbqState
 * @see ParserAbqCLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <algorithm>
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
 * Registers standard Abaqus concentrated-load history records.
 *
 * Each data line uses `node-or-nset, dof, magnitude`. The logical identity is
 * the original target token together with the generalized DOF, so a later
 * `OP=MOD` line replaces the corresponding active definition instead of adding
 * a second FEMaster load to it. An empty `OP=NEW` block is valid and removes all
 * active concentrated loads without defining replacements.
 *
 * Follower and imaginary harmonic loads remain unsupported. `OP` must be
 * consistent across all `*CLOAD` cards in one step.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining active load definitions.
 */
inline void register_cload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("CLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define or modify active Abaqus concentrated nodal loads.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("MOD preserves active CLOADs; NEW replaces the complete CLOAD category")
                .flag("FOLLOWER").doc("Unsupported follower concentrated-load flag")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("CLOAD must appear after a supported procedure inside STEP");
            }
            if (state.procedure == "EIGENFREQ") {
                throw std::runtime_error("CLOAD is not supported in a FREQUENCY step");
            }
            if (keys.has("FOLLOWER")) {
                throw std::runtime_error("CLOAD FOLLOWER is not supported");
            }
            if (keys.has("REAL") && keys.has("IMAGINARY")) {
                throw std::runtime_error("CLOAD REAL and IMAGINARY are mutually exclusive");
            }
            if (keys.has("IMAGINARY")) {
                throw std::runtime_error("CLOAD IMAGINARY is not supported by the real-load harmonic solver");
            }

            const std::string op = keys.raw("OP");
            if (!state.cload_op.empty() && state.cload_op != op) {
                throw std::runtime_error("All CLOAD cards in one STEP must use the same OP value");
            }
            if (state.cload_op.empty()) {
                state.cload_op = op;
                if (op == "NEW") {
                    state.cloads.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            if (!amplitude->empty() && !parser.model()._data->amplitudes.has(*amplitude)) {
                throw std::runtime_error("CLOAD references unknown amplitude '" + *amplitude + "'");
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or node identifier")
                    .one<int>().name("DOF").desc("Abaqus generalized degree of freedom 1--6")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Total load magnitude")
                )
                .bind([&parser, amplitude](const std::string& target,
                                           int dof,
                                           fem::Precision magnitude) {
                    if (dof < 1 || dof > 6) {
                        throw std::runtime_error("CLOAD supports only structural DOFs 1 through 6");
                    }

                    auto& state = parser.abaqus_state();
                    auto record = std::find_if(
                        state.cloads.begin(), state.cloads.end(),
                        [&](const ParserAbqCLoad& item) {
                            return item.target == target && item.dof == dof;
                        }
                    );

                    if (record == state.cloads.end()) {
                        state.cloads.push_back(ParserAbqCLoad{
                            target,
                            dof,
                            magnitude,
                            *amplitude,
                            state.step_index
                        });
                    } else {
                        record->magnitude     = magnitude;
                        record->amplitude     = *amplitude;
                        record->modified_step = state.step_index;
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
