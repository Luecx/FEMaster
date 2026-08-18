/**
 * @file register_boundary.inl
 * @brief Registers Abaqus *BOUNDARY displacement and rotation constraints.
 *
 * Boundary conditions are translated directly into the FEMaster support
 * collector used by the single supported Abaqus analysis step. Definitions at
 * the root become constant model-level constraints and may prescribe nonzero
 * values. Definitions inside `*STEP` may prescribe supported nonzero static
 * values. Node-set targets are expanded while parsing so nodal `*TRANSFORM`
 * coordinate systems can be assigned to each generated support.
 *
 * @see ParserAbqState
 * @see model::Model::add_support
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <charconv>
#include <limits>
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
 * Registers the supported Abaqus displacement boundary-condition syntax.
 *
 * Data use `target, first_dof, last_dof, magnitude`. Omitted `last_dof` selects
 * only the first DOF and an omitted magnitude prescribes zero. Root-level
 * definitions become constant constraints and may prescribe nonzero values
 * without an amplitude. Nonzero values inside the analysis step are supported
 * only by static FEMaster procedures. Named amplitudes for nonzero step-local
 * values are evaluated at the end of a linear static step because FEMaster
 * support equations do not carry amplitudes.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active step and nodal transforms.
 */
inline void register_boundary(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "STEP"}));
        command.doc("Define Abaqus displacement or rotation constraints.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("DISPLACEMENT").allowed({"DISPLACEMENT"})
                    .doc("Only displacement-type boundary conditions are supported")
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(!state.step_active || parser.active_loadcase(),
                "BOUNDARY inside STEP must appear after a supported procedure");

            parser.model()._data->supp_cols.activate("__ABQ_STEP_SUPPORTS");

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "BOUNDARY references unknown amplitude '", *amplitude, "'");
            logging::error(state.step_active || amplitude->empty(),
                "BOUNDARY outside STEP does not support AMPLITUDE");
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

                    logging::error(first_dof >= 1 && first_dof <= 6
                                && last_dof >= first_dof && last_dof <= 6,
                        "BOUNDARY supports only structural DOFs 1 through 6");

                    // Root-level values are constant model constraints; step-local
                    // nonzero values additionally depend on the active procedure.
                    const bool in_step = parser.abaqus_state().step_active;
                    if (in_step && magnitude != fem::Precision(0)) {
                        const std::string& procedure = parser.active_loadcase_type();
                        logging::error(procedure == "LINEARSTATIC"
                                    || procedure == "NONLINEARSTATIC"
                                    || procedure == "STATIC_RIKS",
                            "Nonzero prescribed BOUNDARY values are supported only for static FEMaster procedures");

                        if (!amplitude->empty()) {
                            logging::error(procedure == "LINEARSTATIC",
                                "Nonzero BOUNDARY AMPLITUDE is unsupported because FEMaster constraints are time-independent");
                            magnitude *= parser.model()._data->amplitudes.get(*amplitude)->evaluate(
                                parser.abaqus_state().step_period
                            );
                        }
                    }

                    fem::StaticVector<6> values;
                    values.setConstant(std::numeric_limits<fem::Precision>::quiet_NaN());
                    for (int dof = first_dof; dof <= last_dof; ++dof) {
                        values[dof - 1] = magnitude;
                    }

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    const auto add_to_node = [&](fem::ID node_id) {
                        std::string orientation;
                        const auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            orientation = transform->second;
                        }
                        model.add_support(node_id, values, orientation);
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
                        "BOUNDARY target '", target, "' is not a node set or node id");
                    add_to_node(node_id);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
