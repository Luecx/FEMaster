/**
 * @file register_loadcase_newmark.inl
 * @brief Registers Newmark-beta integration parameters for transient analyses.
 *
 * The `NEWMARK` child command reads beta and gamma for the implicit fixed-step
 * time integrator and stores them on the active `Transient` load case. Default
 * values remain beta = 0.25 and gamma = 0.5 when the command is absent.
 *
 * Time stepping, effective operator construction and state advancement remain
 * responsibilities of the transient solver.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"

#include "../../../loadcase/linear_transient.h"

namespace fem::io::reader::commands {

inline void register_loadcase_newmark(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("NEWMARK", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set Newmark-β integration parameters (β, γ). Defaults are 0.25, 0.5.");

        command.variant(
            fem::io::dsl::Variant::make()
                .doc("One data line: β, γ")
                .segment(
                    fem::io::dsl::Segment::make()
                        .range(fem::io::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::io::dsl::Pattern::make()
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

                            logging::error(false, "NEWMARK not supported for loadcase type " + base->type_name());
                        })
                )
        );
    });
}

} // namespace fem::io::reader::commands
