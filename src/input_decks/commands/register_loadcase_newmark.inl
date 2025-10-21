#pragma once
/**
 * @file register_loadcase_newmark.inl
 * @brief Register *NEWMARK inside *LOADCASE (Transient): (β, γ) parameters.
 */

#include <array>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../core/logging.h"
#include "../../core/types_num.h"

#include "../../loadcase/linear_transient.h"

namespace fem::input_decks::commands {

inline void register_loadcase_newmark(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("NEWMARK", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set Newmark-β integration parameters (β, γ). Defaults are 0.25, 0.5.");

        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: β, γ")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>()
                                .name("NEWMARK")
                                .desc("β, γ parameters for Newmark-β.")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&parser](const std::array<fem::Precision, 2>& bg) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "NEWMARK must appear inside *LOADCASE.");

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                lc->set_newmark(bg[0], bg[1]);
                                return;
                            }

                            logging::error(false, "NEWMARK not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
