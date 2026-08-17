/**
 * @file register_step.inl
 * @brief Registers supported Abaqus analysis steps and procedure cards.
 *
 * Every Abaqus `*STEP` is translated to one independent FEMaster load case.
 * Supported procedures are general/static Riks analysis, eigenfrequency
 * extraction, linear eigenvalue buckling, implicit linear transient dynamics and
 * direct steady-state harmonic dynamics. `*END STEP` attaches the generated
 * step-local load/support collectors, executes the load case and clears the
 * active parser state.
 *
 * FEMaster currently does not propagate converged state between load cases, so
 * Abaqus history continuation and `OP=MOD` persistence across multiple steps are
 * deliberately not emulated here. Each supported step starts from the FEMaster
 * model reference state.
 *
 * @see ParserAbqState
 * @see loadcase::LinearStatic
 * @see loadcase::NonlinearStatic
 * @see loadcase::LinearEigenfrequency
 * @see loadcase::LinearBuckling
 * @see loadcase::Transient
 * @see loadcase::LinearHarmonic
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the Abaqus step block, supported procedure cards and `*END STEP`.
 *
 * `NLGEOM=YES` selects FEMaster nonlinear static analysis for a general static
 * step. `*STATIC, RIKS` selects the arc-length controller regardless of the
 * explicit NLGEOM setting. Procedure data are reduced to the controls that have
 * direct FEMaster counterparts; unsupported termination criteria or adaptive
 * dynamic time stepping are rejected or intentionally not represented rather
 * than introducing new solver semantics in the input reader.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser owning active load-case and step state.
 */
inline void register_step(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    // ---------------------------------------------------------------------
    // STEP block
    // ---------------------------------------------------------------------
    registry.command("STEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin an Abaqus analysis step mapped to one FEMaster load case.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional().doc("Optional Abaqus step name")
                .key("NLGEOM").optional("NO").allowed({"YES", "NO"})
                    .doc("Enable geometric nonlinearity for supported static analysis")
                .key("INC").optional("100").doc("Maximum number of increments")
                .key("AMPLITUDE").optional().allowed({"RAMP", "STEP"})
                    .doc("Default Abaqus step amplitude mode")
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

            state.step_active = true;
            state.step_index  = state.next_step_index++;
            state.max_increments = max_increments;
            state.nlgeom = keys.raw("NLGEOM") == "YES";
            state.step_period = fem::Precision(1);
            state.step_name = keys.has("NAME") ? keys.raw("NAME") : std::string{};
            state.step_amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            state.procedure.clear();
            state.load_collector = "__ABQ_STEP_" + std::to_string(state.step_index) + "_LOADS";
            state.support_collector = "__ABQ_STEP_" + std::to_string(state.step_index) + "_SUPPORTS";
            state.load_collector_used = false;
            state.support_collector_used = false;
        });

        // Reaching the end of the STEP scope while it is still active means the
        // explicit Abaqus END STEP card was missing.
        command.on_exit([&parser](const fem::io::dsl::Keys&) {
            if (parser.abaqus_state().step_active) {
                throw std::runtime_error("Abaqus STEP is missing *END STEP");
            }
        });

        command.variant(fem::io::dsl::Variant::make());
    });

    // ---------------------------------------------------------------------
    // General static and Riks static procedure
    // ---------------------------------------------------------------------
    registry.command("STATIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus general static or Riks analysis to FEMaster static analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("RIKS").doc("Select arc-length path following")
                .flag("DIRECT").doc("Request fixed Abaqus increments")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }

            const bool riks = keys.has("RIKS");
            const int id = parser.next_loadcase_id();

            if (riks || state.nlgeom) {
                auto loadcase = std::make_unique<loadcase::NonlinearStatic>(
                    id, &parser.writer(), &parser.model()
                );

                // Abaqus defaults correspond to one complete normalized step
                // with a minimum automatic increment of 1e-5 of the step scale.
                // Set them here because an omitted STATIC data line has no
                // segment callback in the DSL engine.
                loadcase->max_increments    = state.max_increments;
                loadcase->initial_increment = fem::Precision(1);
                loadcase->minimum_increment = fem::Precision(1e-5);
                loadcase->maximum_increment = fem::Precision(1);
                loadcase->control = riks
                    ? loadcase::NonlinearControl::ArcLength
                    : loadcase::NonlinearControl::LoadControl;

                parser.set_active_loadcase(
                    std::move(loadcase),
                    riks ? "STATIC_RIKS" : "NONLINEARSTATIC"
                );
                state.procedure = riks ? "STATIC_RIKS" : "NONLINEARSTATIC";
            } else {
                parser.set_active_loadcase(
                    std::make_unique<loadcase::LinearStatic>(id, &parser.writer(), &parser.model()),
                    "LINEARSTATIC"
                );
                state.procedure = "LINEARSTATIC";
            }
        });

        // Riks accepts up to eight values, while FEMaster consumes the arc-length
        // increment controls represented by the first four.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_present("RIKS"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 8>().name("DATA")
                        .desc("Riks arc-length controls")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<fem::Precision, 8>& data) {
                    auto* loadcase = parser.active_loadcase_as<loadcase::NonlinearStatic>();
                    if (!loadcase) {
                        throw std::runtime_error("STATIC, RIKS did not create a nonlinear load case");
                    }
                    if (!std::isnan(data[4]) || !std::isnan(data[5]) ||
                        !std::isnan(data[6]) || !std::isnan(data[7])) {
                        throw std::runtime_error(
                            "STATIC, RIKS termination by LPF/displacement is not supported"
                        );
                    }

                    const fem::Precision period = std::isnan(data[1]) || data[1] == fem::Precision(0)
                        ? fem::Precision(1) : data[1];
                    const fem::Precision initial = std::isnan(data[0]) || data[0] == fem::Precision(0)
                        ? period : data[0];
                    const fem::Precision minimum = std::isnan(data[2]) || data[2] == fem::Precision(0)
                        ? std::min(initial, fem::Precision(1e-5) * period) : data[2];
                    const fem::Precision maximum = std::isnan(data[3]) || data[3] == fem::Precision(0)
                        ? period : data[3];

                    // FEMaster's arc-length controls are normalized to an
                    // equivalent load-factor scale, while Abaqus supplies arc
                    // lengths relative to its total arc-length scale.
                    parser.abaqus_state().step_period = period;
                    loadcase->initial_increment = initial / period;
                    loadcase->minimum_increment = minimum / period;
                    loadcase->maximum_increment = maximum / period;
                })
            )
        );

        // General static data map directly to FEMaster load-control increment
        // limits when NLGEOM selects the nonlinear solver. Linear static retains
        // only the step period for final amplitude evaluation.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::key_present("RIKS")
            ))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 4>().name("DATA")
                        .desc("Static time-increment controls")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<fem::Precision, 4>& data) {
                    const fem::Precision period = std::isnan(data[1]) || data[1] == fem::Precision(0)
                        ? fem::Precision(1) : data[1];
                    const fem::Precision initial = std::isnan(data[0]) || data[0] == fem::Precision(0)
                        ? period : data[0];
                    const fem::Precision minimum = std::isnan(data[2]) || data[2] == fem::Precision(0)
                        ? std::min(initial, fem::Precision(1e-5) * period) : data[2];
                    const fem::Precision maximum = std::isnan(data[3]) || data[3] == fem::Precision(0)
                        ? period : data[3];

                    parser.abaqus_state().step_period = period;
                    if (auto* loadcase = parser.active_loadcase_as<loadcase::NonlinearStatic>()) {
                        loadcase->initial_increment = initial / period;
                        loadcase->minimum_increment = minimum / period;
                        loadcase->maximum_increment = maximum / period;
                    }
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Natural-frequency extraction
    // ---------------------------------------------------------------------
    registry.command("FREQUENCY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus frequency extraction to FEMaster eigenfrequency analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("EIGENSOLVER").optional("LANCZOS").allowed({"LANCZOS", "SUBSPACE"})
                    .doc("Accepted Abaqus eigensolver label; FEMaster selects its own backend")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }
            parser.set_active_loadcase(
                std::make_unique<loadcase::LinearEigenfrequency>(
                    parser.next_loadcase_id(), &parser.writer(), &parser.model(), 10
                ),
                "EIGENFREQ"
            );
            state.procedure = "EIGENFREQ";
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<int>().name("NUM").desc("Number of eigenvalues")
                    .fixed<fem::Precision, 5>().name("REST").desc("Abaqus spectral controls ignored by FEMaster")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser](int count, const std::array<fem::Precision, 5>&) {
                    if (count <= 0) {
                        throw std::runtime_error("FREQUENCY requires a positive eigenvalue count");
                    }
                    parser.active_loadcase_as<loadcase::LinearEigenfrequency>()->num_eigenvalues = count;
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Linear eigenvalue buckling
    // ---------------------------------------------------------------------
    registry.command("BUCKLE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus eigenvalue buckling to FEMaster linear buckling analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("EIGENSOLVER").optional("SUBSPACE").allowed({"LANCZOS", "SUBSPACE"})
                    .doc("Accepted Abaqus eigensolver label; FEMaster selects its own backend")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }
            parser.set_active_loadcase(
                std::make_unique<loadcase::LinearBuckling>(
                    parser.next_loadcase_id(), &parser.writer(), &parser.model(), 10
                ),
                "LINEARBUCKLING"
            );
            state.procedure = "LINEARBUCKLING";
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<int>().name("NUM").desc("Number of buckling eigenvalues")
                    .fixed<fem::Precision, 4>().name("REST").desc("Abaqus buckling spectral controls ignored by FEMaster")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser](int count, const std::array<fem::Precision, 4>&) {
                    if (count <= 0) {
                        throw std::runtime_error("BUCKLE requires a positive eigenvalue count");
                    }
                    parser.active_loadcase_as<loadcase::LinearBuckling>()->num_eigenvalues = count;
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Implicit transient dynamics
    // ---------------------------------------------------------------------
    registry.command("DYNAMIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus implicit transient dynamics to FEMaster linear Newmark analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("DIRECT").doc("Abaqus fixed-increment request")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }
            if (state.nlgeom) {
                throw std::runtime_error("DYNAMIC with NLGEOM=YES is not supported by FEMaster's linear transient solver");
            }

            parser.set_active_loadcase(
                std::make_unique<loadcase::Transient>(
                    parser.next_loadcase_id(), &parser.writer(), &parser.model()
                ),
                "LINEARTRANSIENT"
            );
            state.procedure = "LINEARTRANSIENT";
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 4>().name("DATA")
                        .desc("Initial increment, period, minimum increment, maximum increment")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<fem::Precision, 4>& data) {
                    if (std::isnan(data[0]) || data[0] <= fem::Precision(0) ||
                        std::isnan(data[1]) || data[1] <= fem::Precision(0)) {
                        throw std::runtime_error("DYNAMIC requires positive initial increment and step period");
                    }

                    auto* loadcase = parser.active_loadcase_as<loadcase::Transient>();
                    loadcase->dt      = data[0];
                    loadcase->t_start = 0.0;
                    loadcase->t_end   = data[1];
                    parser.abaqus_state().step_period = data[1];
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Direct steady-state harmonic dynamics
    // ---------------------------------------------------------------------
    registry.command("STEADYSTATEDYNAMICS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map direct Abaqus steady-state dynamics to FEMaster harmonic response.");

        auto frequency_scale = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("DIRECT").doc("Require direct physical-DOF harmonic analysis")
                .key("INTERVAL").optional("RANGE").allowed({"RANGE"})
                    .doc("Only direct frequency ranges are supported")
                .key("FREQUENCYSCALE").optional("LOGARITHMIC").allowed({"LOGARITHMIC", "LINEAR"})
                    .doc("Frequency spacing, logarithmic by Abaqus default")
        );

        command.on_enter([&parser, frequency_scale](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }
            if (!keys.has("DIRECT")) {
                throw std::runtime_error("Only STEADY STATE DYNAMICS, DIRECT is supported");
            }
            if (state.nlgeom) {
                throw std::runtime_error("STEADY STATE DYNAMICS does not support NLGEOM in FEMaster");
            }

            *frequency_scale = keys.raw("FREQUENCYSCALE");
            parser.set_active_loadcase(
                std::make_unique<loadcase::LinearHarmonic>(
                    parser.next_loadcase_id(), &parser.writer(), &parser.model()
                ),
                "LINEARHARMONIC"
            );
            state.procedure = "LINEARHARMONIC";
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 5>().name("DATA")
                        .desc("Lower frequency, upper frequency, point count, bias, scale factor")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser, frequency_scale](const std::array<fem::Precision, 5>& data) {
                    auto* loadcase = parser.active_loadcase_as<loadcase::LinearHarmonic>();
                    if (!loadcase || std::isnan(data[0]) || data[0] < fem::Precision(0)) {
                        throw std::runtime_error("STEADY STATE DYNAMICS requires a non-negative lower frequency");
                    }
                    if (!std::isnan(data[3]) && std::abs(data[3] - fem::Precision(1)) > fem::Precision(1e-12)) {
                        throw std::runtime_error("Biased steady-state frequency spacing is not supported");
                    }
                    if (!std::isnan(data[4]) && std::abs(data[4] - fem::Precision(1)) > fem::Precision(1e-12)) {
                        throw std::runtime_error("Steady-state frequency scale factors other than 1 are not supported");
                    }

                    const fem::Precision lower = data[0];
                    const fem::Precision upper = std::isnan(data[1]) ? fem::Precision(0) : data[1];
                    if (upper == fem::Precision(0)) {
                        loadcase->frequencies.push_back(lower);
                        return;
                    }
                    if (upper < lower || std::isnan(data[2])) {
                        throw std::runtime_error("Steady-state frequency range requires upper >= lower and a point count");
                    }

                    const int count = static_cast<int>(std::llround(data[2]));
                    if (count < 2 || std::abs(data[2] - static_cast<fem::Precision>(count)) > fem::Precision(1e-12)) {
                        throw std::runtime_error("Steady-state frequency point count must be an integer >= 2");
                    }

                    if (*frequency_scale == "LOGARITHMIC" && lower <= fem::Precision(0)) {
                        throw std::runtime_error("Logarithmic steady-state frequency ranges require a positive lower frequency");
                    }

                    for (int i = 0; i < count; ++i) {
                        const fem::Precision xi = static_cast<fem::Precision>(i) /
                                                  static_cast<fem::Precision>(count - 1);
                        if (*frequency_scale == "LINEAR") {
                            loadcase->frequencies.push_back(lower + xi * (upper - lower));
                        } else {
                            loadcase->frequencies.push_back(std::exp(
                                std::log(lower) + xi * (std::log(upper) - std::log(lower))
                            ));
                        }
                    }
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Explicit step terminator: attach step-local collectors and execute
    // ---------------------------------------------------------------------
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Finish and execute the currently open Abaqus step.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            auto* base  = parser.active_loadcase();
            if (!state.step_active || !base || state.procedure.empty()) {
                throw std::runtime_error("END STEP requires one active supported procedure");
            }

            // Attach only collectors defined in this Abaqus step. FEMaster load
            // cases are independent, so previous Abaqus history is not inherited.
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
                throw std::runtime_error(std::string("Abaqus STEP execution failed: ") + error.what());
            }

            parser.clear_active_loadcase();
            state.step_active = false;
            state.procedure.clear();
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
