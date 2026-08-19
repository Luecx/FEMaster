/**
 * @file register_step.inl
 * @brief Registers the supported Abaqus analysis step and procedure syntax.
 *
 * The Abaqus reader accepts at most one analysis `*STEP`. Supported procedure
 * cards map directly to FEMaster load-case classes, and step-local loads and
 * boundary conditions are written to dedicated FEMaster collectors while their
 * input cards are parsed.
 *
 * Supported procedures include static load control, Riks arc-length analysis,
 * linear static perturbation, eigenfrequency extraction, linear buckling, direct
 * implicit transient dynamics and direct steady-state harmonic response. Abaqus
 * node/element file and print requests are accepted and ignored because result
 * output remains controlled by FEMaster's writers.
 *
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
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the single supported Abaqus analysis step and its procedures.
 *
 * `*STEP` stores the active step controls and creates dedicated load and support
 * collectors. Exactly one supported procedure card must appear before step-local
 * loads and boundary conditions. Each procedure attaches those collectors to its
 * FEMaster load case, and `*END STEP` executes that load case.
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
        command.doc("Begin the single supported Abaqus analysis step.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").optional().doc("Optional Abaqus step name")
                .key("NLGEOM").optional().allowed({"", "YES", "NO", "OFF"})
                    .doc("Enable geometric nonlinearity when present without a value or set to YES")
                .key("INC").optional("100").doc("Maximum number of increments")
                .key("AMPLITUDE").optional().allowed({"RAMP", "STEP"})
                    .doc("Default amplitude mode for loads defined in this step")
                .flag("PERTURBATION").doc("Linear perturbation procedure")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(!state.step_seen,
                "The FEMaster Abaqus reader supports at most one analysis STEP");
            logging::error(!state.step_active && !parser.active_loadcase(),
                "Nested or unfinished Abaqus STEP blocks are not supported");

            const int max_increments = keys.get<int>("INC");
            logging::error(max_increments > 0,
                "STEP requires INC > 0");

            // Store the normalized step controls. A bare NLGEOM key and
            // NLGEOM=YES enable geometric nonlinearity; NO and OFF disable it.
            state.step_seen      = true;
            state.step_active    = true;
            state.max_increments = max_increments;
            state.nlgeom         = keys.has("NLGEOM") && keys.get<bool>("NLGEOM");
            state.perturbation   = keys.has("PERTURBATION");
            state.step_period    = Precision(1);
            state.step_amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};

            parser.model()._data->load_cols.activate("__ABQ_STEP_LOADS");
            parser.model()._data->supp_cols.activate("__ABQ_STEP_SUPPORTS");
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

            if (!state.perturbation && (riks || state.nlgeom)) {
                auto loadcase = std::make_unique<loadcase::NonlinearStatic>();

                loadcase->loads.push_back("__ABQ_STEP_LOADS");
                loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
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
                parser.begin_loadcase(std::move(loadcase));
            } else {
                auto loadcase = std::make_unique<loadcase::LinearStatic>();
                loadcase->loads.push_back("__ABQ_STEP_LOADS");
                loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
                parser.begin_loadcase(std::move(loadcase));
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
                    auto* loadcase = dynamic_cast<loadcase::NonlinearStatic*>(parser.active_loadcase());
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
                    if (auto* loadcase = dynamic_cast<loadcase::NonlinearStatic*>(parser.active_loadcase())) {
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
            logging::error(parser.abaqus_state().step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");

            auto loadcase = std::make_unique<loadcase::LinearEigenfrequency>();
            loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
            parser.begin_loadcase(std::move(loadcase));
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
                    auto* loadcase = dynamic_cast<loadcase::LinearEigenfrequency*>(parser.active_loadcase());
                    loadcase->num_eigenvalues = count;
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
            logging::error(parser.abaqus_state().step_active && !parser.active_loadcase(),
                "STEP must contain exactly one supported procedure card");

            auto loadcase = std::make_unique<loadcase::LinearBuckling>();
            loadcase->loads.push_back("__ABQ_STEP_LOADS");
            loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
            parser.begin_loadcase(std::move(loadcase));
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
                    auto* loadcase = dynamic_cast<loadcase::LinearBuckling*>(parser.active_loadcase());
                    loadcase->num_eigenvalues = count;
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
                "DYNAMIC with NLGEOM is not supported by FEMaster's linear transient solver");
            logging::error(!state.perturbation,
                "DYNAMIC is not supported in a PERTURBATION step");

            auto loadcase = std::make_unique<loadcase::Transient>();
            loadcase->loads.push_back("__ABQ_STEP_LOADS");
            loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
            parser.begin_loadcase(std::move(loadcase));
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

                    auto* loadcase = dynamic_cast<loadcase::Transient*>(parser.active_loadcase());
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
            auto loadcase = std::make_unique<loadcase::LinearHarmonic>();
            loadcase->loads.push_back("__ABQ_STEP_LOADS");
            loadcase->supps.push_back("__ABQ_STEP_SUPPORTS");
            parser.begin_loadcase(std::move(loadcase));
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
                    auto* loadcase = dynamic_cast<loadcase::LinearHarmonic*>(parser.active_loadcase());
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
    // Ignored Abaqus output requests
    // ---------------------------------------------------------------------
    const auto register_ignored_output_request = [&](const std::string& name) {
        registry.command(name, [](fem::io::dsl::Command& command) {
            command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "STEP"}));
            command.doc("Accept and ignore an Abaqus result-output request.");

            command.variant(fem::io::dsl::Variant::make()
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(0))
                    .pattern(fem::io::dsl::Pattern::make()
                        .fixed<std::string, 64>().name("OUTPUT").desc("Ignored output variables")
                            .on_missing(std::string{"_"}).on_empty(std::string{"_"})
                    )
                    .bind([](const std::array<std::string, 64>&) {})
                )
            );
        });
    };

    register_ignored_output_request("NODEFILE");
    register_ignored_output_request("ELFILE");
    register_ignored_output_request("NODEPRINT");
    register_ignored_output_request("ELPRINT");

    // ---------------------------------------------------------------------
    // Step terminator
    // ---------------------------------------------------------------------
    registry.command("ENDSTEP", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Execute and finish the active Abaqus analysis step.");

        command.on_enter([&parser](const fem::io::dsl::Keys&) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && parser.active_loadcase() != nullptr,
                "END STEP requires one active supported procedure");

            parser.end_loadcase();
            state.step_active = false;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands_abq
