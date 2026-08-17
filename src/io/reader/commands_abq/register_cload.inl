/**
 * @file register_cload.inl
 * @brief Registers step-local Abaqus *CLOAD concentrated nodal loads.
 *
 * Abaqus supplies one generalized nodal degree of freedom and one magnitude per
 * record. FEMaster stores a six-component `CLoad`, so every record is expanded
 * into a sparse six-vector and added to the current step's generated load
 * collector. Node-set targets are expanded to individual nodes so an Abaqus
 * `*TRANSFORM` can be resolved independently for every target node.
 *
 * Named amplitudes are preserved for linear transient and harmonic analysis.
 * Linear static and buckling analyses use the amplitude value at the end of the
 * step because only the final load state is solved. Custom amplitude paths are
 * rejected for nonlinear static/Riks where FEMaster currently uses one
 * proportional reference load vector.
 *
 * @see ParserAbqState
 * @see bc::CLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../bc/amplitude.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the standard Abaqus concentrated-load format.
 *
 * Each line uses `node-or-nset, dof, magnitude`, with DOFs 1--6 mapped to
 * `[Fx,Fy,Fz,Mx,My,Mz]`. Loads are stored node-by-node so any transform assigned
 * by `*TRANSFORM` is copied to the resulting FEMaster load object rather than
 * changing the model's global nodal coordinate basis.
 *
 * `OP=NEW`, follower concentrated loads and imaginary harmonic loads are
 * rejected because the corresponding Abaqus history/follower/complex-load
 * semantics are not represented by the current FEMaster load object.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser providing step and transform state.
 */
inline void register_cload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("CLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Apply Abaqus nodal force or moment records in the current step.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD"})
                    .doc("Only additive/modifying current-step loads are supported")
                .flag("FOLLOWER").doc("Unsupported follower concentrated-load flag")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("CLOAD must appear after a supported procedure inside STEP");
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

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            parser.model()._data->load_cols.activate(state.load_collector);
            state.load_collector_used = true;
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
                    if (dof < 1 || dof > 6) {
                        throw std::runtime_error("CLOAD supports only structural DOFs 1 through 6");
                    }

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();

                    // Transient and harmonic procedures evaluate the named
                    // amplitude at the current time/frequency. Static and
                    // buckling procedures need only the final step value.
                    std::string stored_amplitude;
                    fem::Precision scale = fem::Precision(1);
                    if (!amplitude->empty()) {
                        if (!model._data->amplitudes.has(*amplitude)) {
                            throw std::runtime_error("CLOAD references unknown amplitude '" + *amplitude + "'");
                        }

                        if (state.procedure == "LINEARTRANSIENT" ||
                            state.procedure == "LINEARHARMONIC") {
                            stored_amplitude = *amplitude;
                        } else if (state.procedure == "LINEARSTATIC" ||
                                   state.procedure == "LINEARBUCKLING") {
                            scale = model._data->amplitudes.get(*amplitude)->evaluate(state.step_period);
                        } else if (state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS") {
                            throw std::runtime_error(
                                "CLOAD AMPLITUDE is not supported for nonlinear static/Riks because FEMaster currently uses proportional load-factor control"
                            );
                        }
                    } else if (state.procedure == "LINEARTRANSIENT" && state.step_amplitude == "RAMP") {
                        const std::string generated = "__ABQ_STEP_" + std::to_string(state.step_index) + "_DEFAULT_AMPLITUDE";
                        if (!model._data->amplitudes.has(generated)) {
                            model.define_amplitude(generated, bc::Interpolation::Linear);
                            model.add_amplitude_sample(generated, fem::Precision(0), fem::Precision(0));
                            model.add_amplitude_sample(generated, state.step_period, fem::Precision(1));
                        }
                        stored_amplitude = generated;
                    } else if ((state.procedure == "NONLINEARSTATIC" || state.procedure == "STATIC_RIKS") &&
                               state.step_amplitude == "STEP") {
                        throw std::runtime_error(
                            "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control"
                        );
                    }

                    fem::Vec6 load = fem::Vec6::Zero();
                    load[dof - 1] = scale * magnitude;

                    const auto add_to_node = [&](fem::ID node_id) {
                        std::string orientation;
                        auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            orientation = transform->second;
                        }
                        model.add_cload(node_id, load, orientation, stored_amplitude);
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
                        throw std::runtime_error("CLOAD target '" + target + "' is not a node set or node id");
                    }
                    add_to_node(static_cast<fem::ID>(value));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
