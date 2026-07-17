// register_loadcase_nonlinear.inl — registers NONLINEAR controls within *LOADCASE

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/nonlinear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_nonlinear(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("NONLINEAR", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Configure NONLINEARSTATIC increment and iteration controls.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("INCREMENTS").optional()
                    .doc("Legacy increment count; sets INITIAL_INCREMENT to 1 / INCREMENTS.")
                .key("MAX_INCREMENTS").optional()
                    .doc("Maximum number of accepted nonlinear increments.")
                .key("INITIAL_INCREMENT").optional()
                    .doc("Initial nonlinear increment size.")
                .key("MINIMUM_INCREMENT").optional()
                    .doc("Minimum nonlinear increment size.")
                .key("MAXIMUM_INCREMENT").optional()
                    .doc("Maximum nonlinear increment size.")
                .key("CONTROL").optional("LOAD")
                    .allowed({"LOAD", "ARC_LENGTH"})
                    .doc("Nonlinear path control: LOAD or ARC_LENGTH.")
                .key("ARC_LENGTH_PSI").optional("1.0")
                    .doc("Weighting factor for the load part of the arc-length constraint.")
                .key("ADAPTIVE").optional("ON")
                .key("GROWTH_FACTOR").optional("1.5")
                    .doc("Growth factor for quickly converged nonlinear increments.")
                .key("CUTBACK_FACTOR").optional("0.5")
                    .doc("Reduction factor for rejected or slowly converged increments.")
                .key("FAST_ITERATIONS").optional("6")
                    .doc("Accepted increments converging within this count are enlarged.")
                .key("SLOW_ITERATIONS").optional("10")
                    .doc("Accepted increments requiring at least this count are reduced.")
                .key("MAXIMUM_CUTBACKS").optional("20")
                    .doc("Maximum consecutive cutbacks allowed for one increment.")
                .key("MAXITER").optional("20")
                .key("TOL").optional("1e-8")
                .key("REGULARIZE_ZERO_ROWS").optional("ON")
                .key("REGULARIZATION_ALPHA").optional("1e-4")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("NONLINEAR must appear inside *LOADCASE");
            }

            auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base);
            if (!lc) {
                throw std::runtime_error("NONLINEAR is only supported for NONLINEARSTATIC loadcases");
            }

            if (keys.has("INCREMENTS")) {
                const int increments = keys.get<int>("INCREMENTS");
                if (increments <= 0) {
                    throw std::runtime_error("NONLINEAR requires INCREMENTS > 0");
                }
                lc->initial_increment =
                    Precision(1) / static_cast<Precision>(increments);
            }

            if (keys.has("MAX_INCREMENTS")) {
                lc->max_increments = keys.get<int>("MAX_INCREMENTS");
            }
            if (keys.has("INITIAL_INCREMENT")) {
                lc->initial_increment = keys.get<Precision>("INITIAL_INCREMENT");
            }
            if (keys.has("MINIMUM_INCREMENT")) {
                lc->minimum_increment = keys.get<Precision>("MINIMUM_INCREMENT");
            }
            if (keys.has("MAXIMUM_INCREMENT")) {
                lc->maximum_increment = keys.get<Precision>("MAXIMUM_INCREMENT");
            }

            lc->control = keys.equals("CONTROL", "ARC_LENGTH")
                ? loadcase::NonlinearControl::ArcLength
                : loadcase::NonlinearControl::LoadControl;
            lc->arc_length_psi                       = keys.get<Precision>("ARC_LENGTH_PSI");
            lc->adaptive_increments                  = keys.get<bool>("ADAPTIVE");
            lc->growth_factor                        = keys.get<Precision>("GROWTH_FACTOR");
            lc->cutback_factor                       = keys.get<Precision>("CUTBACK_FACTOR");
            lc->fast_iterations                      = keys.get<int>("FAST_ITERATIONS");
            lc->slow_iterations                      = keys.get<int>("SLOW_ITERATIONS");
            lc->maximum_cutbacks                     = keys.get<int>("MAXIMUM_CUTBACKS");
            lc->max_iterations                       = keys.get<int>("MAXITER");
            lc->tolerance                            = keys.get<Precision>("TOL");
            lc->regularize_zero_stiffness_rows       = keys.get<bool>("REGULARIZE_ZERO_ROWS");
            lc->zero_stiffness_regularization_alpha = keys.get<Precision>("REGULARIZATION_ALPHA");
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
