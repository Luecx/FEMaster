#pragma once
/**
 * @file register_loadcase_write_every.inl
 * @brief Register *WRITE EVERY inside *LOADCASE (Transient): output cadence.
 */

#include <array>
#include <string>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../core/logging.h"
#include "../../core/types_num.h"

#include "../../loadcase/linear_transient.h"

namespace fem::input_decks::commands {

inline void register_loadcase_write_every(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("WRITEEVERY", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Control output frequency during transient analysis. "
            "Use TYPE=STEPS with an integer N (write every N steps), or "
            "TYPE=TIME with a positive Δt_write in seconds (overrides steps)."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").optional("STEPS").allowed({"STEPS","TIME"})
        );

        auto type = std::make_shared<std::string>("STEPS");
        command.on_enter([type](const fem::dsl::Keys& keys) {
            *type = keys.raw("TYPE");
        });

        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: either integer steps (TYPE=STEPS) or Δt_write seconds (TYPE=TIME).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .one<fem::Precision>()
                                .name("VALUE")
                                .desc("N (steps) or Δt_write (seconds), depending on TYPE.")
                        )
                        .bind([&parser, type](fem::Precision value) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "WRITE EVERY must appear inside *LOADCASE.");

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                if (*type == "TIME") {
                                    lc->set_write_every_time(static_cast<double>(value));
                                } else {
                                    // Treat as step count; clamp to >=1
                                    const int steps = std::max<int>(1, static_cast<int>(std::llround(value)));
                                    lc->set_write_every(steps);
                                }
                                return;
                            }

                            logging::error(false, "WRITE EVERY not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
