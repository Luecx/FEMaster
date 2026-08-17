/**
 * @file register_history.inl
 * @brief Materializes persistent Abaqus load/BC definitions at *END STEP.
 *
 * FEMaster load cases remain mechanically independent, while Abaqus load and
 * boundary-condition definitions are propagated logically between steps by the
 * dedicated command handlers. This adapter converts the complete active snapshot
 * into fresh FEMaster load/support collectors immediately before the current
 * load case is executed.
 *
 * Step-time amplitudes on records modified in the current step are active only
 * for that step. After successful execution their end value is stored as the
 * constant propagated magnitude used by later steps unless the record is
 * modified again. For an explicitly requested default `RAMP`, a transient load
 * redefinition is decomposed into its previous constant total value plus a
 * linearly ramped difference to the new total value.
 *
 * @see ParserAbqState
 * @see commands_abq::register_step
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "../parser_abq.h"
#include "../../dsl/registry.h"
#include "../../../bc/amplitude.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Replaces the structural `*END STEP` callback from `register_step` with logical
 * Abaqus history materialization and FEMaster load-case execution.
 *
 * Only definition history is persistent. Previous FEMaster collectors are never
 * attached to a later load case; instead, every step receives a new collector
 * containing the current active snapshot after `OP=MOD`/`OP=NEW` processing.
 * Nodal targets are expanded only here so `*TRANSFORM` can be resolved per node.
 *
 * @param registry Final-stage DSL registry containing the ENDSTEP command.
 * @param parser Abaqus parser owning logical history and active load-case state.
 */
inline void register_history(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.doc("Materialize active Abaqus load/BC definitions and execute the independent FEMaster step.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& model = parser.model();
            auto& state = parser.abaqus_state();
            auto* base  = parser.active_loadcase();
            if (!state.step_active || !base || state.procedure.empty()) {
                throw std::runtime_error("END STEP requires one active supported procedure");
            }

            const bool nonlinear = state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS";
            const bool transient = state.procedure == "LINEARTRANSIENT";
            const bool use_loads = state.procedure != "EIGENFREQ";

            // Explicit named amplitudes can be retained only by procedures that
            // evaluate loads repeatedly. Static/buckling calculations need only
            // the end-of-step total value; nonlinear proportional loading cannot
            // reproduce an arbitrary independent amplitude history.
            const auto resolve_explicit_amplitude = [&](const std::string& amplitude,
                                                         int modified_step)
                -> std::pair<Precision, std::string> {
                if (modified_step != state.step_index || amplitude.empty()) {
                    return {Precision(1), std::string{}};
                }
                if (!model._data->amplitudes.has(amplitude)) {
                    throw std::runtime_error("Unknown Abaqus amplitude '" + amplitude + "'");
                }

                if (state.procedure == "LINEARTRANSIENT" ||
                    state.procedure == "LINEARHARMONIC") {
                    return {Precision(1), amplitude};
                }
                if (state.procedure == "LINEARSTATIC" ||
                    state.procedure == "LINEARBUCKLING") {
                    return {
                        model._data->amplitudes.get(amplitude)->evaluate(state.step_period),
                        std::string{}
                    };
                }
                if (nonlinear) {
                    throw std::runtime_error(
                        "Named load AMPLITUDE is not supported for nonlinear static/Riks proportional loading"
                    );
                }
                return {Precision(1), std::string{}};
            };

            // Lazily create the one unit ramp shared by every load that is
            // redefined without an explicit amplitude in this transient step.
            const auto ramp_amplitude = [&]() -> std::string {
                const std::string name =
                    "__ABQ_STEP_" + std::to_string(state.step_index) + "_DEFAULT_AMPLITUDE";
                if (!model._data->amplitudes.has(name)) {
                    model.define_amplitude(name, bc::Interpolation::Linear);
                    model.add_amplitude_sample(name, Precision(0), Precision(0));
                    model.add_amplitude_sample(name, state.step_period, Precision(1));
                }
                return name;
            };

            // -----------------------------------------------------------------
            // Boundary-condition snapshot
            // -----------------------------------------------------------------
            // Always attach an explicit step-local collector, even when it is
            // empty. FEMaster interprets an empty support-name list as the global
            // ALL collector, which would otherwise reactivate supports created
            // for earlier Abaqus steps after BOUNDARY, OP=NEW.
            model._data->supp_cols.activate(state.support_collector);
            state.support_collector_used = true;

            for (const auto& record : state.boundaries) {
                Precision magnitude = record.magnitude;

                if (record.modified_step == state.step_index &&
                    !record.amplitude.empty() && magnitude != Precision(0)) {
                    if (!model._data->amplitudes.has(record.amplitude)) {
                        throw std::runtime_error(
                            "BOUNDARY references unknown amplitude '" + record.amplitude + "'"
                        );
                    }
                    if (state.procedure != "LINEARSTATIC") {
                        throw std::runtime_error(
                            "Nonzero BOUNDARY AMPLITUDE is unsupported because FEMaster constraints are time-independent"
                        );
                    }
                    magnitude *= model._data->amplitudes.get(record.amplitude)->evaluate(state.step_period);
                }

                if (magnitude != Precision(0) &&
                    state.procedure != "LINEARSTATIC" &&
                    state.procedure != "NONLINEARSTATIC" &&
                    state.procedure != "STATIC_RIKS") {
                    throw std::runtime_error(
                        "Nonzero prescribed BOUNDARY values are supported only for static FEMaster procedures"
                    );
                }

                StaticVector<6> values;
                values.setConstant(std::numeric_limits<Precision>::quiet_NaN());
                values[record.dof - 1] = magnitude;

                const auto add_to_node = [&](ID node_id) {
                    std::string orientation;
                    const auto transform = state.node_transforms.find(node_id);
                    if (transform != state.node_transforms.end()) {
                        orientation = transform->second;
                    }
                    model.add_support(node_id, values, orientation);
                };

                if (model._data->node_sets.has(record.target)) {
                    for (const ID node_id : *model._data->node_sets.get(record.target)) {
                        add_to_node(node_id);
                    }
                } else {
                    std::size_t parsed = 0;
                    const long value = std::stol(record.target, &parsed);
                    if (parsed != record.target.size()) {
                        throw std::runtime_error(
                            "BOUNDARY target '" + record.target + "' is not a node set or node id"
                        );
                    }
                    add_to_node(static_cast<ID>(value));
                }
            }

            // -----------------------------------------------------------------
            // Load snapshot
            // -----------------------------------------------------------------
            if (use_loads &&
                (!state.cloads.empty() || !state.dloads.empty() || !state.dsloads.empty())) {
                model._data->load_cols.activate(state.load_collector);
                state.load_collector_used = true;

                for (const auto& record : state.cloads) {
                    const auto add_record = [&](Precision magnitude,
                                                const std::string& amplitude) {
                        if (magnitude == Precision(0)) {
                            return;
                        }

                        Vec6 load = Vec6::Zero();
                        load[record.dof - 1] = magnitude;

                        const auto add_to_node = [&](ID node_id) {
                            std::string orientation;
                            const auto transform = state.node_transforms.find(node_id);
                            if (transform != state.node_transforms.end()) {
                                orientation = transform->second;
                            }
                            model.add_cload(node_id, load, orientation, amplitude);
                        };

                        if (model._data->node_sets.has(record.target)) {
                            for (const ID node_id : *model._data->node_sets.get(record.target)) {
                                add_to_node(node_id);
                            }
                        } else {
                            std::size_t parsed = 0;
                            const long value = std::stol(record.target, &parsed);
                            if (parsed != record.target.size()) {
                                throw std::runtime_error(
                                    "CLOAD target '" + record.target + "' is not a node set or node id"
                                );
                            }
                            add_to_node(static_cast<ID>(value));
                        }
                    };

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && transient &&
                        state.step_amplitude == "RAMP") {
                        add_record(record.previous_magnitude, std::string{});
                        add_record(
                            record.magnitude - record.previous_magnitude,
                            ramp_amplitude()
                        );
                        continue;
                    }

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && nonlinear &&
                        state.step_amplitude == "STEP") {
                        throw std::runtime_error(
                            "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control"
                        );
                    }

                    const auto [scale, amplitude] =
                        resolve_explicit_amplitude(record.amplitude, record.modified_step);
                    add_record(scale * record.magnitude, amplitude);
                }

                for (const auto& record : state.dloads) {
                    const auto add_record = [&](Precision magnitude,
                                                const std::string& amplitude) {
                        if (magnitude == Precision(0)) {
                            return;
                        }

                        const Vec3 gravity = magnitude * Vec3{
                            record.direction[0], record.direction[1], record.direction[2]
                        };

                        if (model._data->elem_sets.has(record.target)) {
                            model.add_vload(record.target, gravity, "", amplitude);
                        } else {
                            std::size_t parsed = 0;
                            const long value = std::stol(record.target, &parsed);
                            if (parsed != record.target.size()) {
                                throw std::runtime_error(
                                    "DLOAD target '" + record.target + "' is not an element set or element id"
                                );
                            }
                            model.add_vload(static_cast<ID>(value), gravity, "", amplitude);
                        }
                    };

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && transient &&
                        state.step_amplitude == "RAMP") {
                        add_record(record.previous_magnitude, std::string{});
                        add_record(
                            record.magnitude - record.previous_magnitude,
                            ramp_amplitude()
                        );
                        continue;
                    }

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && nonlinear &&
                        state.step_amplitude == "STEP") {
                        throw std::runtime_error(
                            "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control"
                        );
                    }

                    const auto [scale, amplitude] =
                        resolve_explicit_amplitude(record.amplitude, record.modified_step);
                    add_record(scale * record.magnitude, amplitude);
                }

                for (const auto& record : state.dsloads) {
                    const auto add_record = [&](Precision magnitude,
                                                const std::string& amplitude) {
                        if (magnitude == Precision(0)) {
                            return;
                        }

                        if (record.type == "P") {
                            if (nonlinear) {
                                throw std::runtime_error(
                                    "DSLOAD P is a follower pressure in Abaqus and is not supported in nonlinear FEMaster steps"
                                );
                            }
                            model.add_pload(record.surface, magnitude, amplitude);
                            return;
                        }

                        if (record.type != "TRVEC") {
                            throw std::runtime_error("Unsupported propagated DSLOAD type '" + record.type + "'");
                        }
                        if (nonlinear && record.follower != "NO") {
                            throw std::runtime_error(
                                "Nonlinear DSLOAD TRVEC requires FOLLOWER=NO; follower traction is not implemented"
                            );
                        }

                        Vec3 direction{
                            record.direction[0], record.direction[1], record.direction[2]
                        };
                        const Precision norm = direction.norm();
                        if (!(norm > Precision(0)) || !std::isfinite(norm)) {
                            throw std::runtime_error("Stored DSLOAD TRVEC direction is invalid");
                        }
                        direction /= norm;

                        model.add_dload(
                            record.surface,
                            magnitude * direction,
                            record.orientation,
                            amplitude
                        );
                    };

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && transient &&
                        state.step_amplitude == "RAMP") {
                        add_record(record.previous_magnitude, std::string{});
                        add_record(
                            record.magnitude - record.previous_magnitude,
                            ramp_amplitude()
                        );
                        continue;
                    }

                    if (record.modified_step == state.step_index &&
                        record.amplitude.empty() && nonlinear &&
                        state.step_amplitude == "STEP") {
                        throw std::runtime_error(
                            "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control"
                        );
                    }

                    const auto [scale, amplitude] =
                        resolve_explicit_amplitude(record.amplitude, record.modified_step);
                    add_record(scale * record.magnitude, amplitude);
                }
            }

            // Attach only the fresh snapshot collector to this independent load
            // case. The explicit support collector is attached even when empty
            // so the native global-ALL fallback is suppressed for Abaqus steps.
            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->supps.push_back(state.support_collector);
            } else if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                lc->supps.push_back(state.support_collector);
            } else if (auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base)) {
                lc->supps.push_back(state.support_collector);
            } else if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                lc->supps.push_back(state.support_collector);
            } else if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                lc->supps.push_back(state.support_collector);
            } else if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                lc->supps.push_back(state.support_collector);
            }

            if (state.load_collector_used) {
                if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                    lc->loads.push_back(state.load_collector);
                } else if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                    lc->loads.push_back(state.load_collector);
                } else if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                    lc->loads.push_back(state.load_collector);
                } else if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                    lc->loads.push_back(state.load_collector);
                } else if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                    lc->loads.push_back(state.load_collector);
                } else {
                    throw std::runtime_error("Loads are not supported in the selected Abaqus procedure");
                }
            }

            try {
                base->run();
            } catch (const std::exception& error) {
                parser.clear_active_loadcase();
                state.step_active = false;
                state.procedure.clear();
                throw std::runtime_error(std::string("Abaqus STEP execution failed: ") + error.what());
            }

            // Freeze current-step explicit amplitudes to the end value used when
            // the same logical definition propagates unchanged into the next
            // step. Default RAMP records already store their new final total.
            Precision end_coordinate = state.step_period;
            if (auto* harmonic = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                if (!harmonic->frequencies.empty()) {
                    end_coordinate = harmonic->frequencies.back();
                }
            }

            const auto freeze_load_records = [&](auto& records) {
                for (auto& record : records) {
                    if (record.modified_step != state.step_index) {
                        continue;
                    }
                    if (!record.amplitude.empty()) {
                        record.magnitude *=
                            model._data->amplitudes.get(record.amplitude)->evaluate(end_coordinate);
                        record.amplitude.clear();
                    }
                    record.previous_magnitude = Precision(0);
                    record.modified_step = 0;
                }
            };

            freeze_load_records(state.cloads);
            freeze_load_records(state.dloads);
            freeze_load_records(state.dsloads);

            for (auto& record : state.boundaries) {
                if (record.modified_step != state.step_index) {
                    continue;
                }
                if (!record.amplitude.empty()) {
                    record.magnitude *=
                        model._data->amplitudes.get(record.amplitude)->evaluate(end_coordinate);
                    record.amplitude.clear();
                }
                record.modified_step = 0;
            }

            parser.clear_active_loadcase();
            state.step_active = false;
            state.procedure.clear();
        });
    });
}

} // namespace fem::io::reader::commands_abq
