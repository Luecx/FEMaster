/**
 * @file register_boundary.inl
 * @brief Registers step-local Abaqus *BOUNDARY displacement constraints.
 *
 * The direct Abaqus boundary format is converted to FEMaster `Support` objects.
 * Node sets are expanded to individual nodes so `*TRANSFORM` can supply the
 * local translational/rotational basis for each constrained node. Only structural
 * displacement/rotation DOFs 1--6 are represented.
 *
 * FEMaster supports fixed prescribed values but does not currently carry a
 * time-dependent amplitude on constraint equations. Named amplitudes are
 * therefore reduced to their end-of-step value for linear static analysis and
 * rejected for nonzero constraints in transient, harmonic and nonlinear path
 * analyses where that reduction would alter the requested history.
 *
 * @see ParserAbqState
 * @see bc::Support
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus direct displacement boundary-condition records.
 *
 * Data use `target, first_dof, last_dof, magnitude`; omitted `last_dof` selects
 * only the first DOF and omitted magnitude prescribes zero. The same magnitude
 * is applied to every DOF in the inclusive range. Abaqus transformed DOFs are
 * represented by assigning the corresponding FEMaster coordinate system to each
 * generated support rather than modifying global nodal coordinates.
 *
 * In Abaqus/Standard a prescribed displacement without an explicit amplitude is
 * ramped over the step. FEMaster nonlinear static analysis applies prescribed
 * support values proportionally with the load factor, so that default maps
 * directly. `OP=NEW` remains unsupported because cross-step history removal is
 * not represented by the current independent load-case mapping.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser providing step and transform state.
 */
inline void register_boundary(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Apply Abaqus displacement/rotation boundary conditions in the current step.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("DISPLACEMENT").allowed({"DISPLACEMENT"})
                    .doc("Only displacement-type boundary conditions are supported")
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD"})
                    .doc("Only additive/modifying current-step constraints are supported")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("BOUNDARY must appear after a supported procedure inside STEP");
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            parser.model()._data->supp_cols.activate(state.support_collector);
            state.support_collector_used = true;
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

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    fem::Precision scale = fem::Precision(1);

                    if (!amplitude->empty() && magnitude != fem::Precision(0)) {
                        if (!model._data->amplitudes.has(*amplitude)) {
                            throw std::runtime_error("BOUNDARY references unknown amplitude '" + *amplitude + "'");
                        }

                        if (state.procedure == "LINEARSTATIC") {
                            scale = model._data->amplitudes.get(*amplitude)->evaluate(state.step_period);
                        } else {
                            throw std::runtime_error(
                                "Nonzero BOUNDARY AMPLITUDE is unsupported because FEMaster constraints are currently time-independent"
                            );
                        }
                    }

                    if (magnitude != fem::Precision(0) &&
                        state.procedure != "LINEARSTATIC" &&
                        state.procedure != "NONLINEARSTATIC" &&
                        state.procedure != "STATIC_RIKS") {
                        throw std::runtime_error(
                            "Nonzero prescribed BOUNDARY values are supported only for static FEMaster procedures"
                        );
                    }

                    fem::StaticVector<6> values;
                    values.setConstant(std::numeric_limits<fem::Precision>::quiet_NaN());
                    for (int dof = first_dof; dof <= last_dof; ++dof) {
                        values[dof - 1] = scale * magnitude;
                    }

                    const auto add_to_node = [&](fem::ID node_id) {
                        std::string orientation;
                        auto transform = state.node_transforms.find(node_id);
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

                    std::size_t parsed = 0;
                    const long value = std::stol(target, &parsed);
                    if (parsed != target.size()) {
                        throw std::runtime_error("BOUNDARY target '" + target + "' is not a node set or node id");
                    }
                    add_to_node(static_cast<fem::ID>(value));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
