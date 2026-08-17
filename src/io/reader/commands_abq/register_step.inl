/**
 * @file register_step.inl
 * @brief Registers supported Abaqus analysis steps and procedure cards.
 *
 * Every supported `*STEP` maps to one mechanically independent FEMaster load
 * case. No displacement, stress, material, contact or other converged solver
 * state is inherited from the preceding step. Load and boundary-condition
 * definition history is handled separately by `register_history`.
 *
 * Supported procedures are static/general and Riks analysis, independent linear
 * static perturbation, eigenfrequency extraction, linear eigenvalue buckling,
 * direct implicit linear transient dynamics and direct steady-state harmonic
 * dynamics.
 *
 * @see commands_abq::register_history
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
 * Registers Abaqus step boundaries and supported procedure definitions.
 *
 * `NLGEOM` applies only to the current FEMaster load case. `*STATIC, RIKS`
 * selects arc-length control, while a normal static step with `NLGEOM=YES`
 * selects nonlinear load control. `*STEP, PERTURBATION` with `*STATIC` maps to
 * an independent `LinearStatic` calculation rather than a perturbation about a
 * previous solved state.
 *
 * Abaqus automatic implicit-dynamic incrementation has no direct FEMaster
 * counterpart, so only `*DYNAMIC, DIRECT` is accepted. Direct steady-state
 * dynamics maps to FEMaster's direct harmonic sweep.
 *
 * `*END STEP` is registered here structurally; its execution callback is supplied
 * by `register_history` in the final parser pass after logical load/BC records
 * have been updated.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser owning the active load-case state.
 */
inline void register_step(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    // ---------------------------------------------------------------------
    // STEP block
    // ---------------------------------------------------------------------
    registry.command("STEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin one mechanically independent FEMaster load case.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional().doc("Optional Abaqus step name")
                .key("NLGEOM").optional("NO").allowed({"YES", "NO"})
                    .doc("Enable geometric nonlinearity for this step only")
                .key("INC").optional("100").doc("Maximum number of increments")
                .key("AMPLITUDE").optional().allowed({"RAMP", "STEP"})
                    .doc("Default amplitude mode for loads modified in this step")
                .flag("PERTURBATION").doc("Independent linear perturbation procedure")
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

            state.step_active     = true;
            state.step_index      = state.next_step_index++;
            state.max_increments  = max_increments;
            state.nlgeom          = keys.raw("NLGEOM") == "YES";
            state.perturbation    = keys.has("PERTURBATION");
            state.step_period     = Precision(1);
            state.step_name       = keys.has("NAME") ? keys.raw("NAME") : std::string{};
            state.step_amplitude  = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            state.procedure.clear();

            state.cload_op.clear();
            state.boundary_op.clear();
            state.dload_op.clear();
            state.dsload_op.clear();

            state.load_collector    = "__ABQ_STEP_" + std::to_string(state.step_index) + "_LOADS";
            state.support_collector = "__ABQ_STEP_" + std::to_string(state.step_index) + "_SUPPORTS";
            state.load_collector_used    = false;
            state.support_collector_used = false;
        });

        command.on_exit([&parser](const fem::io::dsl::Keys&) {
            if (parser.abaqus_state().step_active) {
                throw std::runtime_error("Abaqus STEP is missing *END STEP");
            }
        });

        command.variant(fem::io::dsl::Variant::make());
    });

    // ---------------------------------------------------------------------
    // Static/general, independent perturbation and Riks
    // ---------------------------------------------------------------------
    registry.command("STATIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus static/general, perturbation or Riks analysis to FEMaster static analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("RIKS").doc("Select arc-length path following")
                .flag("DIRECT").doc("Use a fixed nonlinear load increment")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }

            const bool riks   = keys.has("RIKS");
            const bool direct = keys.has("DIRECT");
            if (riks && state.perturbation) {
                throw std::runtime_error("STATIC, RIKS cannot be used in a PERTURBATION step");
            }

            const int id = parser.next_loadcase_id();
            if (!state.perturbation && (riks || state.nlgeom)) {
                auto loadcase = std::make_unique<loadcase::NonlinearStatic>(
                    id, &parser.writer(), &parser.model()
                );

                loadcase->max_increments      = state.max_increments;
                loadcase->initial_increment   = Precision(1);
                loadcase->minimum_increment   = Precision(1e-5);
                loadcase->maximum_increment   = riks
                    ? std::numeric_limits<Precision>::max()
                    : Precision(1);
                loadcase->adaptive_increments = !direct;
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

        // Riks: first four entries are the arc-length increment controls. The
        // remaining Abaqus LPF/displacement termination criteria are not mapped.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_present("RIKS"))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<Precision, 8>().name("DATA")
                        .desc("Riks arc-length controls")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<Precision, 8>& data) {
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

                    const Precision period = std::isnan(data[1]) || data[1] == Precision(0)
                        ? Precision(1) : data[1];
                    const Precision initial = std::isnan(data[0]) || data[0] == Precision(0)
                        ? period : data[0];
                    const Precision minimum = std::isnan(data[2]) || data[2] == Precision(0)
                        ? std::min(initial, Precision(1e-5) * period) : data[2];
                    const bool maximum_omitted = std::isnan(data[3]) || data[3] == Precision(0);

                    parser.abaqus_state().step_period = period;
                    loadcase->initial_increment = initial / period;
                    loadcase->minimum_increment = minimum / period;
                    loadcase->maximum_increment = maximum_omitted
                        ? std::numeric_limits<Precision>::max()
                        : data[3] / period;
                })
            )
        );

        // General static: increment controls matter only for nonlinear FEMaster
        // load control. Linear static retains the period for amplitude evaluation.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::negate(
                fem::io::dsl::Condition::key_present("RIKS")
            ))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<Precision, 4>().name("DATA")
                        .desc("Static time-increment controls")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<Precision, 4>& data) {
                    if (parser.abaqus_state().perturbation) {
                        for (const Precision value : data) {
                            if (!std::isnan(value)) {
                                throw std::runtime_error(
                                    "STATIC in a PERTURBATION step does not accept general-step increment data"
                                );
                            }
                        }
                        return;
                    }

                    const Precision period = std::isnan(data[1]) || data[1] == Precision(0)
                        ? Precision(1) : data[1];
                    const Precision initial = std::isnan(data[0]) || data[0] == Precision(0)
                        ? period : data[0];
                    const Precision minimum = std::isnan(data[2]) || data[2] == Precision(0)
                        ? std::min(initial, Precision(1e-5) * period) : data[2];
                    const Precision maximum = std::isnan(data[3]) || data[3] == Precision(0)
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
                    .doc("Accepted label; FEMaster selects its own backend")
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
                    .fixed<Precision, 5>().name("REST").desc("Unused Abaqus spectral controls")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser](int count, const std::array<Precision, 5>&) {
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
        command.doc("Map Abaqus eigenvalue buckling to independent FEMaster linear buckling.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("EIGENSOLVER").optional("SUBSPACE").allowed({"LANCZOS", "SUBSPACE"})
                    .doc("Accepted label; FEMaster selects its own backend")
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
                    .fixed<Precision, 4>().name("REST").desc("Unused Abaqus buckling controls")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser](int count, const std::array<Precision, 4>&) {
                    if (count <= 0) {
                        throw std::runtime_error("BUCKLE requires a positive eigenvalue count");
                    }
                    parser.active_loadcase_as<loadcase::LinearBuckling>()->num_eigenvalues = count;
                })
            )
        );
    });

    // ---------------------------------------------------------------------
    // Direct implicit transient dynamics
    // ---------------------------------------------------------------------
    registry.command("DYNAMIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus DYNAMIC, DIRECT to FEMaster linear Newmark analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("DIRECT").doc("Required fixed-increment Abaqus dynamics mode")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || parser.active_loadcase()) {
                throw std::runtime_error("STEP must contain exactly one supported procedure card");
            }
            if (!keys.has("DIRECT")) {
                throw std::runtime_error(
                    "Only DYNAMIC, DIRECT is supported because FEMaster transient analysis uses a fixed time increment"
                );
            }
            if (state.nlgeom) {
                throw std::runtime_error("DYNAMIC with NLGEOM=YES is not supported by FEMaster's linear transient solver");
            }
            if (state.perturbation) {
                throw std::runtime_error("DYNAMIC is not supported in an independent PERTURBATION step");
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
                    .fixed<Precision, 4>().name("DATA")
                        .desc("Fixed increment, period, unused minimum, unused maximum")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser](const std::array<Precision, 4>& data) {
                    if (std::isnan(data[0]) || data[0] <= Precision(0) ||
                        std::isnan(data[1]) || data[1] <= Precision(0)) {
                        throw std::runtime_error("DYNAMIC, DIRECT requires positive increment and step period");
                    }

                    auto* loadcase = parser.active_loadcase_as<loadcase::Transient>();
                    loadcase->dt      = data[0];
                    loadcase->t_start = Precision(0);
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
                    .fixed<Precision, 5>().name("DATA")
                        .desc("Lower frequency, upper frequency, point count, bias, scale factor")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser, frequency_scale](const std::array<Precision, 5>& data) {
                    auto* loadcase = parser.active_loadcase_as<loadcase::LinearHarmonic>();
                    if (!loadcase || std::isnan(data[0]) || data[0] < Precision(0)) {
                        throw std::runtime_error("STEADY STATE DYNAMICS requires a non-negative lower frequency");
                    }

                    const Precision bias = std::isnan(data[3]) ? Precision(1) : data[3];
                    if (std::abs(bias - Precision(1)) > Precision(1e-12)) {
                        throw std::runtime_error("Biased steady-state frequency spacing is not supported");
                    }
                    if (!std::isnan(data[4]) && std::abs(data[4] - Precision(1)) > Precision(1e-12)) {
                        throw std::runtime_error("Steady-state frequency scale factors other than 1 are not supported");
                    }

                    const Precision lower = data[0];
                    const Precision upper = std::isnan(data[1]) ? Precision(0) : data[1];
                    if (upper == Precision(0)) {
                        loadcase->frequencies.push_back(lower);
                        return;
                    }
                    if (upper < lower) {
                        throw std::runtime_error("Steady-state frequency range requires upper >= lower");
                    }

                    int count = 20;
                    if (!std::isnan(data[2])) {
                        const int requested = static_cast<int>(std::llround(data[2]));
                        if (std::abs(data[2] - static_cast<Precision>(requested)) > Precision(1e-12)) {
                            throw std::runtime_error("Steady-state frequency point count must be an integer");
                        }
                        if (requested >= 2) {
                            count = requested;
                        }
                    }

                    if (*frequency_scale == "LOGARITHMIC" && lower <= Precision(0)) {
                        throw std::runtime_error("Logarithmic steady-state frequency ranges require a positive lower frequency");
                    }

                    for (int i = 0; i < count; ++i) {
                        const Precision xi = static_cast<Precision>(i) /
                                             static_cast<Precision>(count - 1);
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

    // ENDSTEP is executed by register_history in the final parser pass.
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Finish the active Abaqus step.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
