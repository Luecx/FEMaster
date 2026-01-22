#pragma once
/**
 * @file register_loadcase_time.inl
 * @brief Register *TIME inside *LOADCASE (Transient): fixed step and end time.
 */

#include <array>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../core/logging.h"
#include "../../core/types_num.h"

#include "../../loadcase/linear_transient.h"

namespace fem::input_decks::commands {

inline void register_loadcase_time(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("TIME", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set time window and step for transient analysis. Accepts either: (t_start, t_end, dt) or (t_end, dt).");

        // Variant A: three values → t_start, t_end, dt
        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: t_start, t_end, dt")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 3>()
                                .name("TIME")
                                .desc("t_start, t_end, dt")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&parser](const std::array<fem::Precision, 3>& T) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "TIME must appear inside *LOADCASE.");

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                const auto t_start = T[0];
                                const auto t_end   = T[1];
                                const auto dt      = T[2];
                                lc->set_time(dt, t_start, t_end);
                                return;
                            }

                            logging::error(false, "TIME not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );

        // Variant B: two values → t_end, dt (t_start defaults to 0)
        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: t_end, dt")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>()
                                .name("TIME")
                                .desc("t_end, dt")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&parser](const std::array<fem::Precision, 2>& T) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "TIME must appear inside *LOADCASE.");

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                const auto t_end = T[0];
                                const auto dt    = T[1];
                                lc->set_time(dt, 0.0, t_end);
                                return;
                            }

                            logging::error(false, "TIME not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
