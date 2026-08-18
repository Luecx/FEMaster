/**
 * @file register_history.inl
 * @brief Materializes Abaqus load and boundary-condition history at *END STEP.
 *
 * The registration converts the active logical Abaqus CLOAD, DLOAD, DSLOAD and
 * BOUNDARY records into FEMaster load and support collectors for the current
 * analysis step. It resolves nodal transforms, load amplitudes, default transient
 * ramps and procedure-specific restrictions before executing the active load case.
 *
 * After execution, step-time amplitude definitions are reduced to their final
 * values for propagation into later steps. The propagated records describe load
 * and boundary-condition definitions only; FEMaster load cases retain independent
 * mechanical initial states.
 *
 * @see ParserAbqState
 * @see commands_abq::register_step
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <charconv>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <system_error>
#include <utility>

#include "../parser_abq.h"
#include "../../dsl/registry.h"
#include "../../../bc/amplitude.h"
#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus `*END STEP` finalization for the supported history syntax.
 *
 * The active logical load and boundary definitions are materialized into the
 * collector names assigned to the open step. Node-set targets are expanded to
 * individual nodes when a nodal `*TRANSFORM` must be assigned. Explicit
 * amplitudes and the default transient `RAMP` are mapped to FEMaster amplitude
 * objects where the selected procedure supports time-dependent loading.
 *
 * @param registry Final-stage DSL registry containing the ENDSTEP command.
 * @param parser Abaqus parser owning logical history and active load-case state.
 */
inline void register_history(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.doc("Materialize active Abaqus load/BC definitions and execute the current FEMaster step.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& model = parser.model();
            auto& state = parser.abaqus_state();
            auto* base  = parser.active_loadcase();

            logging::error(state.step_active && base != nullptr && !state.procedure.empty(),
                "END STEP requires one active supported procedure");

            const bool nonlinear = state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS";
            const bool transient = state.procedure == "LINEARTRANSIENT";
            const bool use_loads = state.procedure != "EIGENFREQ";

            const auto parse_id = [](const std::string& value, const char* label) -> ID {
                ID id{};
                const char* begin = value.data();
                const char* end   = begin + value.size();
                const auto [ptr, ec] = std::from_chars(begin, end, id);
                logging::error(ec == std::errc{} && ptr == end,
                    label, " '", value, "' is not a valid numeric identifier");
                return id;
            };

            // Resolve explicit amplitudes for the selected analysis procedure.
            const auto resolve_explicit_amplitude = [&](const std::string& amplitude,
                                                         int modified_step)
                -> std::pair<Precision, std::string> {
                if (modified_step != state.step_index || amplitude.empty()) {
                    return {Precision(1), std::string{}};
                }

                logging::error(model._data->amplitudes.has(amplitude),
                    "Unknown Abaqus amplitude '", amplitude, "'");

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

                logging::error(!nonlinear,
                    "Named load AMPLITUDE is not supported for nonlinear static/Riks proportional loading");
                return {Precision(1), std::string{}};
            };

            // Unit ramp used for transient definitions without an explicit amplitude.
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
            if (!state.boundaries.empty()) {
                model._data->supp_cols.activate(state.support_collector);
                state.support_collector_used = true;

                for (const auto& record : state.boundaries) {
                    Precision magnitude = record.magnitude;

                    if (record.modified_step == state.step_index &&
                        !record.amplitude.empty() && magnitude != Precision(0)) {
                        logging::error(model._data->amplitudes.has(record.amplitude),
                            "BOUNDARY references unknown amplitude '", record.amplitude, "'");
                        logging::error(state.procedure == "LINEARSTATIC",
                            "Nonzero BOUNDARY AMPLITUDE is unsupported because FEMaster constraints are time-independent");
                        magnitude *= model._data->amplitudes.get(record.amplitude)->evaluate(state.step_period);
                    }

                    logging::error(magnitude == Precision(0)
                                || state.procedure == "LINEARSTATIC"
                                || state.procedure == "NONLINEARSTATIC"
                                || state.procedure == "STATIC_RIKS",
                        "Nonzero prescribed BOUNDARY values are supported only for static FEMaster procedures");

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
                        add_to_node(parse_id(record.target, "BOUNDARY target"));
                    }
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
                            add_to_node(parse_id(record.target, "CLOAD target"));
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

                    logging::error(!(record.modified_step == state.step_index
                                  && record.amplitude.empty() && nonlinear
                                  && state.step_amplitude == "STEP"),
                        "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control");

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
                            model.add_vload(
                                parse_id(record.target, "DLOAD target"),
                                gravity,
                                "",
                                amplitude
                            );
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

                    logging::error(!(record.modified_step == state.step_index
                                  && record.amplitude.empty() && nonlinear
                                  && state.step_amplitude == "STEP"),
                        "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control");

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
                            logging::error(!nonlinear,
                                "DSLOAD P is a follower pressure in Abaqus and is not supported in nonlinear FEMaster steps");
                            model.add_pload(record.surface, magnitude, amplitude);
                            return;
                        }

                        logging::error(record.type == "TRVEC",
                            "Unsupported propagated DSLOAD type '", record.type, "'");
                        logging::error(!nonlinear || record.follower == "NO",
                            "Nonlinear DSLOAD TRVEC requires FOLLOWER=NO; follower traction is not implemented");

                        Vec3 direction{
                            record.direction[0], record.direction[1], record.direction[2]
                        };
                        const Precision norm = direction.norm();
                        logging::error(norm > Precision(0) && std::isfinite(norm),
                            "Stored DSLOAD TRVEC direction is invalid");
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

                    logging::error(!(record.modified_step == state.step_index
                                  && record.amplitude.empty() && nonlinear
                                  && state.step_amplitude == "STEP"),
                        "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control");

                    const auto [scale, amplitude] =
                        resolve_explicit_amplitude(record.amplitude, record.modified_step);
                    add_record(scale * record.magnitude, amplitude);
                }
            }

            // -----------------------------------------------------------------
            // Collector assignment
            // -----------------------------------------------------------------
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
                    logging::error(false,
                        "Loads are not supported in the selected Abaqus procedure");
                }
            }

            // -----------------------------------------------------------------
            // Execute load case
            // -----------------------------------------------------------------
            try {
                base->run();
            } catch (const std::exception& error) {
                parser.clear_active_loadcase();
                state.step_active = false;
                state.procedure.clear();
                logging::error(false,
                    "Abaqus STEP execution failed: ", error.what());
            }

            // -----------------------------------------------------------------
            // Prepare propagated definitions for the next step
            // -----------------------------------------------------------------
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
