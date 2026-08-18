/**
 * @file register_step.inl
 * @brief Registers supported Abaqus analysis-step and procedure syntax.
 *
 * The registration maps Abaqus `*STEP` blocks and supported procedure cards to
 * FEMaster load-case classes. Supported procedures include static load control,
 * Riks arc-length analysis, linear static perturbation, eigenfrequency extraction,
 * linear buckling, direct implicit transient dynamics and direct steady-state
 * harmonic response.
 *
 * Each step defines one FEMaster load case. Step parameters such as `NLGEOM`,
 * increment controls and frequency ranges are transferred to the corresponding
 * load-case object. Load and boundary-condition definitions are handled by the
 * dedicated history registrations.
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
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus step boundaries and the supported analysis procedures.
 *
 * `*STEP` stores step-local control parameters. `*STATIC` selects linear static,
 * nonlinear load-control or Riks analysis according to `PERTURBATION`, `NLGEOM`
 * and `RIKS`. `*FREQUENCY`, `*BUCKLE`, `*DYNAMIC, DIRECT` and
 * `*STEADY STATE DYNAMICS, DIRECT` create the corresponding FEMaster load cases.
 * `*END STEP` is registered as the structural terminator and its execution is
 * supplied by `register_history` in the analysis-data stage.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser owning the active step and load-case state.
 */
inline void register_step(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    // ---------------------------------------------------------------------
    // STEP block
    // ---------------------------------------------------------------------
    registry.command("STEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin one FEMaster load case from an Abaqus STEP block.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional().doc("Optional Abaqus step name")
                .key("NLGEOM").optional("NO").allowed({"YES", "NO"})
                    .doc("Enable geometric nonlinearity for this step")
                .key("INC").optional("100").doc("Maximum number of increments")
                .key("AMPLITUDE").optional().allowed({"RAMP", "STEP"})
                    .doc("Default amplitude mode for loads defined in this step")
                .flag("PERTURBATION").doc("Linear perturbation procedure")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(!state.step_active && !parser.active_loadcase(),
                "Nested or unfinished Abaqus STEP blocks are not supported");

            const int max_increments = keys.get<int>("INC");
            logging::error(max_increments > 0,
                "STEP requires INC > 0");

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
            logging::error(!parser.abaqus_state().step_active,
                "Abaqus STEP is missing *END STEP");
        });

        command.variant(fem::io::dsl::Variant::make());
    });

    // ---------------------------------------------------------------------
    // Static, perturbation and Riks procedures
    // ---------------------------------------------------------------------
    registry.command("STATIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Map Abaqus static, perturbation or Riks analysis to FEMaster static analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("RIKS").doc("Select arc-length path following")
                .flag("DIRECT").doc("Use a fixed nonlinear load increment")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");

            const bool riks   = keys.has("RIKS");
            const bool direct = keys.has("DIRECT");
            logging::error(!(riks && state.perturbation),
                "STATIC, RIKS cannot be used in a PERTURBATION step");

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

        // Riks arc-length increment controls
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
                    logging::error(loadcase != nullptr,
                        "STATIC, RIKS did not create a nonlinear load case");
                    logging::error(std::isnan(data[4]) && std::isnan(data[5])
                                && std::isnan(data[6]) && std::isnan(data[7]),
                        "STATIC, RIKS termination by LPF/displacement is not supported");

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

        // Static time and increment controls
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
                            logging::error(std::isnan(value),
                                "STATIC in a PERTURBATION step does not accept general-step increment data");
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
            logging::error(state.step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");

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
                    logging::error(count > 0,
                        "FREQUENCY requires a positive eigenvalue count");
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
                    .doc("Accepted label; FEMaster selects its own backend")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");

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
                    logging::error(count > 0,
                        "BUCKLE requires a positive eigenvalue count");
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
            logging::error(state.step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");
            logging::error(keys.has("DIRECT"),
                "Only DYNAMIC, DIRECT is supported because FEMaster transient analysis uses a fixed time increment");
            logging::error(!state.nlgeom,
                "DYNAMIC with NLGEOM=YES is not supported by FEMaster's linear transient solver");
            logging::error(!state.perturbation,
                "DYNAMIC is not supported in a PERTURBATION step");

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
                    logging::error(!std::isnan(data[0]) && data[0] > Precision(0)
                                && !std::isnan(data[1]) && data[1] > Precision(0),
                        "DYNAMIC, DIRECT requires positive increment and step period");

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
            logging::error(state.step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");
            logging::error(keys.has("DIRECT"),
                "Only STEADY STATE DYNAMICS, DIRECT is supported");
            logging::error(!state.nlgeom,
                "STEADY STATE DYNAMICS does not support NLGEOM in FEMaster");

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
                    logging::error(loadcase != nullptr && !std::isnan(data[0]) && data[0] >= Precision(0),
                        "STEADY STATE DYNAMICS requires a non-negative lower frequency");

                    const Precision bias = std::isnan(data[3]) ? Precision(1) : data[3];
                    logging::error(std::abs(bias - Precision(1)) <= Precision(1e-12),
                        "Biased steady-state frequency spacing is not supported");
                    logging::error(std::isnan(data[4])
                                || std::abs(data[4] - Precision(1)) <= Precision(1e-12),
                        "Steady-state frequency scale factors other than 1 are not supported");

                    const Precision lower = data[0];
                    const Precision upper = std::isnan(data[1]) ? Precision(0) : data[1];
                    if (upper == Precision(0)) {
                        loadcase->frequencies.push_back(lower);
                        return;
                    }

                    logging::error(upper >= lower,
                        "Steady-state frequency range requires upper >= lower");

                    int count = 20;
                    if (!std::isnan(data[2])) {
                        const int requested = static_cast<int>(std::llround(data[2]));
                        logging::error(std::abs(data[2] - static_cast<Precision>(requested)) <= Precision(1e-12),
                            "Steady-state frequency point count must be an integer");
                        if (requested >= 2) {
                            count = requested;
                        }
                    }

                    logging::error(*frequency_scale != "LOGARITHMIC" || lower > Precision(0),
                        "Logarithmic steady-state frequency ranges require a positive lower frequency");

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

    // ---------------------------------------------------------------------
    // Step terminator
    // ---------------------------------------------------------------------
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Finish the active Abaqus step.");
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
