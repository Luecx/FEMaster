/**
 * @file register_boundary.inl
 * @brief Registers step-dependent Abaqus *BOUNDARY definitions.
 *
 * Boundary conditions are retained as logical target/DOF records across steps.
 * `OP=MOD` updates or adds those records; `OP=NEW` clears the complete active
 * boundary-condition set before the current step supplies its replacements.
 *
 * The original target spelling is preserved until `*END STEP`, where node sets
 * are expanded and nodal `*TRANSFORM` coordinate systems are copied to the
 * generated FEMaster supports.
 *
 * @see ParserAbqState
 * @see ParserAbqBoundary
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
 * Registers direct Abaqus displacement boundary-condition history records.
 *
 * Data use `target, first_dof, last_dof, magnitude`. Ranges are split into one
 * logical record per DOF, making subsequent `OP=MOD` replacement unambiguous.
 * Time-dependent prescribed displacements are validated when the active snapshot
 * is materialized because FEMaster constraints do not currently carry amplitude
 * objects.
 *
 * `OP` must be consistent across all `*BOUNDARY` cards in one step.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining active boundary definitions.
 */
inline void register_boundary(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define or modify active Abaqus displacement boundary conditions.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("DISPLACEMENT").allowed({"DISPLACEMENT"})
                    .doc("Only displacement-type boundary conditions are supported")
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("MOD preserves active boundaries; NEW replaces the complete boundary set")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("BOUNDARY must appear after a supported procedure inside STEP");
            }

            const std::string op = keys.raw("OP");
            if (!state.boundary_op.empty() && state.boundary_op != op) {
                throw std::runtime_error("All BOUNDARY cards in one STEP must use the same OP value");
            }
            if (state.boundary_op.empty()) {
                state.boundary_op = op;
                if (op == "NEW") {
                    state.boundaries.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            if (!amplitude->empty() && !parser.model()._data->amplitudes.has(*amplitude)) {
                throw std::runtime_error("BOUNDARY references unknown amplitude '" + *amplitude + "'");
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or node identifier")
                    .one<int>().name("FIRST_DOF").desc("First constrained DOF")
                    .one<int>().name("LAST_DOF").desc("Last constrained DOF")
                        .on_missing(-1).on_empty(-1)
                    .one<fem::Precision>().name("MAGNITUDE").desc("Prescribed value")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&parser, amplitude](const std::string& target,
                                           int first_dof,
                                           int last_dof,
                                           fem::Precision magnitude) {
                    if (last_dof < 0) {
                        last_dof = first_dof;
                    }
                    if (first_dof < 1 || first_dof > 6 ||
                        last_dof < first_dof || last_dof > 6) {
                        throw std::runtime_error("BOUNDARY supports only structural DOFs 1 through 6");
                    }

                    auto& state = parser.abaqus_state();
                    for (int dof = first_dof; dof <= last_dof; ++dof) {
                        auto record = std::find_if(
                            state.boundaries.begin(), state.boundaries.end(),
                            [&](const ParserAbqBoundary& item) {
                                return item.target == target && item.dof == dof;
                            }
                        );

                        if (record == state.boundaries.end()) {
                            state.boundaries.push_back(ParserAbqBoundary{
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
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
