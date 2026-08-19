/**
 * @file register_loadcase_frequency.inl
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

#pragma once

#include <array>
#include <cmath>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../loadcase/linear_harmonic.h"

namespace fem::io::reader::commands {

inline void register_loadcase_frequency(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("FREQUENCIES", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
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

                            auto* lc = dynamic_cast<fem::loadcase::LinearHarmonic*>(base);
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
}

} // namespace fem::io::reader::commands
