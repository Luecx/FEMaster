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
        command.doc("Set fixed time step Δt and end time t_end for the transient analysis.");

        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: Δt, t_end")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>()
                                .name("TIME")
                                .desc("Δt (time step), t_end (total duration).")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&parser](const std::array<fem::Precision, 2>& dt_T) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "TIME must appear inside *LOADCASE.");

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                lc->set_time(dt_T[0], dt_T[1]);
                                return;
                            }

                            logging::error(false, "TIME not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
