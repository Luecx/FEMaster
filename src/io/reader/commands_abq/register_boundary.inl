/**
 * @file register_boundary.inl
 * @brief Registers Abaqus *BOUNDARY displacement and rotation constraints.
 *
 * Boundary conditions are stored as logical Abaqus target/DOF records and are
 * propagated between analysis steps according to `OP=MOD` and `OP=NEW`. Ranged
 * input is split into one record per generalized degree of freedom so individual
 * constraints can be updated deterministically.
 *
 * The original node or node-set target remains available until `*END STEP`, where
 * node sets are expanded and nodal `*TRANSFORM` coordinate systems are assigned
 * to the generated FEMaster supports.
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
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the supported Abaqus displacement boundary-condition syntax.
 *
 * Data use `target, first_dof, last_dof, magnitude`. Omitted `last_dof` selects
 * only the first DOF and an omitted magnitude prescribes zero. `OP=MOD` preserves
 * and updates the active boundary set, while `OP=NEW` clears the complete set
 * before new records are read.
 *
 * Time-dependent nonzero prescribed values are validated during snapshot
 * materialization because FEMaster support equations do not carry amplitudes.
 * All BOUNDARY cards within one step must use the same `OP` value.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active boundary-definition state.
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
            logging::error(state.step_active && parser.active_loadcase(),
                "BOUNDARY must appear after a supported procedure inside STEP");

            const std::string op = keys.raw("OP");
            logging::error(state.boundary_op.empty() || state.boundary_op == op,
                "All BOUNDARY cards in one STEP must use the same OP value");
            if (state.boundary_op.empty()) {
                state.boundary_op = op;
                if (op == "NEW") {
                    state.boundaries.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "BOUNDARY references unknown amplitude '", *amplitude, "'");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
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

                    logging::error(first_dof >= 1 && first_dof <= 6
                                && last_dof >= first_dof && last_dof <= 6,
                        "BOUNDARY supports only structural DOFs 1 through 6");

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
