#pragma once
/**
 * @file register_loadcase_frequency.inl
 * @brief Register *FREQUENCY inside a linear harmonic load case.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../loadcase/linear_harmonic.h"

namespace fem::io::reader::commands {

inline void register_loadcase_frequency(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("FREQUENCY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Define excitation frequencies for a linear harmonic response analysis.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SCALE").optional().allowed({"LINEAR", "LIST"})
        );

        auto scale = std::make_shared<std::string>("LINEAR");
        command.on_enter([scale](const fem::io::dsl::Keys& keys) {
            if (keys.has("SCALE")) {
                *scale = keys.raw("SCALE");
            }
        });

        command.variant(
            fem::io::dsl::Variant::make()
                .doc("LINEAR: one line start, end, number_of_points. LIST: one frequency per line.")
                .segment(
                    fem::io::dsl::Segment::make()
                        .range(fem::io::dsl::LineRange{}.min(1))
                        .pattern(
                            fem::io::dsl::Pattern::make()
                                .dynamic<fem::Precision>()
                                .name("FREQUENCIES")
                                .desc("Frequency definition in cycles per unit time.")
                        )
                        .bind([&parser, scale](const std::vector<fem::Precision>& values) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "FREQUENCY must appear inside *LOADCASE.");

                            auto* lc = dynamic_cast<fem::loadcase::LinearHarmonic*>(base);
                            logging::error(lc != nullptr,
                                           "FREQUENCY not supported for loadcase type " +
                                           parser.active_loadcase_type());

                            lc->frequencies.clear();

                            if (*scale == "LIST") {
                                logging::error(!values.empty(), "FREQUENCY LIST must not be empty.");
                                lc->frequencies = values;
                                return;
                            }

                            logging::error(values.size() == 3,
                                           "FREQUENCY LINEAR requires start, end, number_of_points.");

                            const fem::Precision start = values[0];
                            const fem::Precision end   = values[1];
                            const int points = static_cast<int>(std::llround(values[2]));

                            logging::error(points >= 1, "FREQUENCY number_of_points must be >= 1.");
                            logging::error(end >= start, "FREQUENCY end must be >= start.");

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
