/**
 * @file register_loadcase_frequency.cpp
 * @brief Registers excitation-frequency sweeps for harmonic response.
 *
 * `FREQUENCIES` defines a linearly spaced sequence from the supplied start and
 * end frequencies and stores it on the active `LinearHarmonic` load case. The
 * command validates the interval and requested point count before constructing
 * the explicit frequency vector consumed by the solver.
 *
 * Dynamic matrix assembly and the response solve at every frequency remain in
 * the harmonic analysis implementation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <cmath>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../loadcase/linear_harmonic.h"

namespace fem::io::reader::commands {

/**
 * Registers excitation-frequency sweeps and optional modal reduction.
 *
 * FREQUENCIES stores an explicit list of non-negative excitation frequencies,
 * linearly spaced between its endpoints. MODALBASIS only enables the Ritz-space
 * solution; eigenvector computation and static-residual enrichment take place
 * within LinearHarmonic::run() after the structural matrices are assembled.
 *
 * @param registry DSL command registry.
 * @param parser Parser providing the active harmonic load case.
 */
void register_loadcase_frequency(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("FREQUENCIES", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Define a linearly spaced excitation-frequency sweep.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SCALE").optional().allowed({"LINEAR"})
        );

        command.variant(
            fem::io::dsl::Variant::make()
                .doc("One data line: start_frequency, end_frequency, number_of_points.")
                .segment(
                    fem::io::dsl::Segment::make()
                        .range(fem::io::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::io::dsl::Pattern::make()
                                .fixed<fem::Precision, 3>()
                                .name("FREQUENCIES")
                                .desc("start, end, number_of_points")
                        )
                        .bind([&parser](const std::array<fem::Precision, 3>& values) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr,
                                "FREQUENCIES must appear inside *LOADCASE.");

                            auto* lc = base->as<fem::loadcase::LinearHarmonic>();
                            logging::error(lc != nullptr,
                                "FREQUENCIES not supported for loadcase type " +
                                base->type_name());

                            const fem::Precision start = values[0];
                            const fem::Precision end   = values[1];
                            const int points = static_cast<int>(std::llround(values[2]));

                            logging::error(std::isfinite(start) && std::isfinite(end),
                                "FREQUENCIES bounds must be finite.");
                            logging::error(start >= 0.0,
                                "FREQUENCIES start must be non-negative.");
                            logging::error(end >= start,
                                "FREQUENCIES end must be >= start.");
                            logging::error(points >= 1,
                                "FREQUENCIES number_of_points must be >= 1.");
                            logging::error(std::abs(values[2] - fem::Precision(points)) < 1e-9,
                                "FREQUENCIES number_of_points must be an integer.");

                            lc->frequencies.resize(points);
                            if (points == 1) {
                                lc->frequencies[0] = start;
                                return;
                            }

                            const fem::Precision increment = (end - start) / fem::Precision(points - 1);
                            for (int i = 0; i < points; ++i) {
                                lc->frequencies[i] = start + fem::Precision(i) * increment;
                            }
                        })
                )
        );
    });

    // A bare MODALBASIS card selects the reduced harmonic solution.
    // The Ritz basis is assembled later from modes and static load responses.
    registry.command("MODALBASIS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Use eigenmodes below twice the maximum excitation frequency and "
            "static load residuals; an empty basis uses the direct solver.");

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            // Validate the active procedure and reject unsupported parameters
            auto* loadcase = parser.active_loadcase();
            logging::error(loadcase != nullptr,
                "MODALBASIS must appear inside *LOADCASE");

            auto* harmonic = loadcase->as<loadcase::LinearHarmonic>();
            logging::error(harmonic != nullptr,
                "MODALBASIS not supported for loadcase type ", loadcase->type_name());
            logging::error(keys._kv.empty(),
                "MODALBASIS does not accept parameters");

            // Only record the choice; construction belongs to the load case
            harmonic->modal_basis = true;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
