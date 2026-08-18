/**
 * @file register_cload.inl
 * @brief Registers Abaqus *CLOAD concentrated nodal load definitions.
 *
 * This registration stores concentrated loads as logical Abaqus records keyed by
 * their original node or node-set target and generalized degree of freedom. The
 * records persist across analysis steps according to `OP=MOD` and `OP=NEW` and
 * are converted to FEMaster `CLoad` objects during `*END STEP` processing.
 *
 * Multiple definitions with the same identity inside one step are stored as
 * separate additive load conditions. Nodal `*TRANSFORM` assignments are resolved
 * when the logical records are materialized for the active load case.
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
#include <string>

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
 * Each data line uses `node-or-nset, dof, magnitude`. `OP=MOD` preserves the
 * active CLOAD set and modifies a uniquely matching propagated definition;
 * `OP=NEW` clears the complete CLOAD set before new records are read. Repeated
 * definitions of the same identity within one step are additive.
 *
 * Follower concentrated loads and imaginary harmonic load components are not
 * supported. All CLOAD cards within one step must use the same `OP` value.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active load-definition state.
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
            logging::error(state.step_active && parser.active_loadcase(),
                "CLOAD must appear after a supported procedure inside STEP");
            logging::error(state.procedure != "EIGENFREQ",
                "CLOAD is not supported in a FREQUENCY step");
            logging::error(!keys.has("FOLLOWER"),
                "CLOAD FOLLOWER is not supported");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "CLOAD REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "CLOAD IMAGINARY is not supported by the real-load harmonic solver");

            const std::string op = keys.raw("OP");
            logging::error(state.cload_op.empty() || state.cload_op == op,
                "All CLOAD cards in one STEP must use the same OP value");
            if (state.cload_op.empty()) {
                state.cload_op = op;
                if (op == "NEW") {
                    state.cloads.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "CLOAD references unknown amplitude '", *amplitude, "'");
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
                    logging::error(dof >= 1 && dof <= 6,
                        "CLOAD supports only structural DOFs 1 through 6");

                    auto& state = parser.abaqus_state();
                    auto first = state.cloads.end();
                    int matches = 0;
                    bool defined_this_step = false;

                    for (auto it = state.cloads.begin(); it != state.cloads.end(); ++it) {
                        if (it->target != target || it->dof != dof) {
                            continue;
                        }
                        if (first == state.cloads.end()) {
                            first = it;
                        }
                        ++matches;
                        defined_this_step = defined_this_step || it->modified_step == state.step_index;
                    }

                    // Repeated definitions in one step represent independent,
                    // additive load conditions with the same Abaqus identity.
                    if (defined_this_step || matches == 0) {
                        state.cloads.push_back(ParserAbqCLoad{
                            target,
                            dof,
                            magnitude,
                            Precision(0),
                            *amplitude,
                            state.step_index
                        });
                        return;
                    }

                    logging::error(matches == 1,
                        "Multiple active CLOADs use target '", target,
                        "' and the same DOF; use OP=NEW to redefine them");

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
