/**
 * @file register_history.inl
 * @brief Applies FEMaster's independent-step semantics to Abaqus load/BC history.
 *
 * Abaqus load and boundary-condition definitions are propagated between steps,
 * while FEMaster mechanical solution state is deliberately not. This adapter
 * overrides only the `*STEP` entry and `*END STEP` finalization hooks installed
 * by `register_step`: each step starts from the FEMaster reference state, but its
 * load/support collectors are materialized from the current logical Abaqus
 * definition snapshot after `OP=MOD`/`OP=NEW` processing.
 *
 * Step-time amplitudes attached to records modified in the current step are used
 * during that step. After successful execution their end value becomes the
 * constant propagated magnitude for later steps, matching Abaqus propagation of
 * unchanged step-time-amplitude loads without carrying solver state.
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
#include "../../dsl/keyword.h"
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
 * Overrides Abaqus step entry/finalization with FEMaster's independent mechanical
 * state and persistent logical load/BC definition history.
 *
 * The function is registered after `register_step` in the final parser pass. It
 * deliberately reuses that command's variants/procedure parsing and replaces
 * only the two hooks whose semantics depend on cross-step history.
 *
 * @param registry Final-stage DSL registry containing the STEP commands.
 * @param parser Abaqus parser owning logical history and active load-case state.
 */
inline void register_history(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    // ---------------------------------------------------------------------
    // Independent STEP entry
    // ---------------------------------------------------------------------
    registry.command("STEP", [&](fem::io::dsl::Command& command) {
        command.doc(
            "Begin one mechanically independent FEMaster load case while preserving Abaqus load/BC definitions."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional().doc("Optional Abaqus step name")
                .key("NLGEOM").optional("NO").allowed({"YES", "NO"})
                    .doc("Enable geometric nonlinearity for this step only")
                .key("INC").optional("100").doc("Maximum number of increments")
                .key("AMPLITUDE").optional().allowed({"RAMP", "STEP"})
                    .doc("Default amplitude mode for loads modified in this step")
                .flag("PERTURBATION").doc("Linear perturbation procedure mapped independently")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("Nested or unfinished Abaqus STEP blocks are not supported");
            }

            const int max_increments = std::stoi(keys.raw("INC"));
            if (max_increments <= 0) {
                throw std::runtime_error("STEP requires INC > 0");
            }

            // Mechanical state is intentionally local to this step. In
            // particular NLGEOM is not inherited from any preceding Abaqus step.
            state.step_active     = true;
            state.step_index      = state.next_step_index++;
            state.max_increments  = max_increments;
            state.nlgeom          = keys.raw("NLGEOM") == "YES";
            state.perturbation    = keys.has("PERTURBATION");
            state.step_period     = Precision(1);
            state.step_name       = keys.has("NAME") ? keys.raw("NAME") : std::string{};
            state.step_amplitude  = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            state.procedure.clear();

            // OP applies independently to each Abaqus load/BC category. The
            // active logical definitions themselves deliberately survive.
            state.cload_op.clear();
            state.boundary_op.clear();
            state.dload_op.clear();
            state.dsload_op.clear();

            state.load_collector    = "__ABQ_STEP_" + std::to_string(state.step_index) + "_LOADS";
            state.support_collector = "__ABQ_STEP_" + std::to_string(state.step_index) + "_SUPPORTS";
            state.load_collector_used    = false;
            state.support_collector_used = false;
        });
    });

    // ---------------------------------------------------------------------
    // Snapshot materialization and execution at END STEP
    // ---------------------------------------------------------------------
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.doc("Materialize the active Abaqus definitions and execute the independent FEMaster step.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& model = parser.model();
            auto& state = parser.abaqus_state();
            auto* base  = parser.active_loadcase();
            if (!state.step_active || !base || state.procedure.empty()) {
                throw std::runtime_error("END STEP requires one active supported procedure");
            }

            const bool nonlinear = state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS";
            const bool use_loads = state.procedure != "EIGENFREQ";

            // A record modified in the current step uses its current amplitude.
            // Propagated records have already been frozen to their prior end
            // value and therefore materialize without an amplitude.
            const auto resolve_load_amplitude = [&](const std::string& amplitude,
                                                     int modified_step)
                -> std::pair<Precision, std::string> {
                if (modified_step != state.step_index) {
                    return {Precision(1), std::string{}};
                }

                if (!amplitude.empty()) {
                    if (!model._data->amplitudes.has(amplitude)) {
                        throw std::runtime_error("Unknown propagated Abaqus amplitude '" + amplitude + "'");
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
                }

                // FEMaster transient analysis can represent the current-step
                // RAMP default directly as one generated linear amplitude.
                if (state.procedure == "LINEARTRANSIENT" && state.step_amplitude == "RAMP") {
                    const std::string generated =
                        "__ABQ_STEP_" + std::to_string(state.step_index) + "_DEFAULT_AMPLITUDE";
                    if (!model._data->amplitudes.has(generated)) {
                        model.define_amplitude(generated, bc::Interpolation::Linear);
                        model.add_amplitude_sample(generated, Precision(0), Precision(0));
                        model.add_amplitude_sample(generated, state.step_period, Precision(1));
                    }
                    return {Precision(1), generated};
                }

                if (nonlinear && state.step_amplitude == "STEP") {
                    throw std::runtime_error(
                        "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control"
                    );
                }

                return {Precision(1), std::string{}};
            };

            // -----------------------------------------------------------------
            // Materialize active boundary conditions
            // -----------------------------------------------------------------
            if (!state.boundaries.empty()) {
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
            }

            // -----------------------------------------------------------------
            // Materialize the complete active load snapshot
            // -----------------------------------------------------------------
            if (use_loads &&
                (!state.cloads.empty() || !state.dloads.empty() || !state.dsloads.empty())) {
                model._data->load_cols.activate(state.load_collector);
                state.load_collector_used = true;

                for (const auto& record : state.cloads) {
                    const auto [scale, amplitude] =
                        resolve_load_amplitude(record.amplitude, record.modified_step);

                    Vec6 load = Vec6::Zero();
                    load[record.dof - 1] = scale * record.magnitude;

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
                }

                for (const auto& record : state.dloads) {
                    const auto [scale, amplitude] =
                        resolve_load_amplitude(record.amplitude, record.modified_step);
                    const Vec3 gravity = scale * record.magnitude * Vec3{
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
                }

                for (const auto& record : state.dsloads) {
                    const auto [scale, amplitude] =
                        resolve_load_amplitude(record.amplitude, record.modified_step);

                    if (record.type == "P") {
                        if (nonlinear) {
                            throw std::runtime_error(
                                "DSLOAD P is a follower pressure in Abaqus and is not supported in nonlinear FEMaster steps"
                            );
                        }
                        model.add_pload(record.surface, scale * record.magnitude, amplitude);
                        continue;
                    }

                    if (record.type != "TRVEC") {
                        throw std::runtime_error("Unsupported propagated DSLOAD type '" + record.type + "'");
                    }
                    if (nonlinear && record.follower != "NO") {
                        throw std::runtime_error(
                            "Nonlinear DSLOAD TRVEC requires FOLLOWER=NO; follower traction is not implemented"
                        );
                    }

                    const Vec3 direction{
                        record.direction[0], record.direction[1], record.direction[2]
                    };
                    model.add_dload(
                        record.surface,
                        scale * record.magnitude * direction,
                        record.orientation,
                        amplitude
                    );
                }
            }

            // Attach exactly the newly materialized snapshot. Previous FEMaster
            // collectors remain stored but are never inherited by this load case.
            if (state.support_collector_used) {
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

            // Step-time amplitude references do not restart when an unchanged
            // load propagates to the next Abaqus step. Freeze every definition
            // modified here to its end-of-step value and clear the amplitude
            // before the next logical snapshot is materialized.
            Precision end_coordinate = state.step_period;
            if (auto* harmonic = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                if (!harmonic->frequencies.empty()) {
                    end_coordinate = harmonic->frequencies.back();
                }
            }

            const auto freeze_records = [&](auto& records) {
                for (auto& record : records) {
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
            };

            freeze_records(state.cloads);
            freeze_records(state.boundaries);
            freeze_records(state.dloads);
            freeze_records(state.dsloads);

            parser.clear_active_loadcase();
            state.step_active = false;
            state.procedure.clear();
        });
    });
}

} // namespace fem::io::reader::commands_abq
