#pragma once
/**
 * @file register_loadcase_damping.inl
 * @brief Register *DAMPING inside *LOADCASE (Transient): Rayleigh C = αM + βK.
 */

#include <array>
#include <string>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../core/logging.h"
#include "../../core/types_num.h"
#include "../../core/types_eig.h"

#include "../../loadcase/linear_transient.h"
#include "../../damping/rayleigh.h"

namespace fem::input_decks::commands {

inline void register_loadcase_damping(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("DAMPING", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Assign damping for the active loadcase. "
            "Currently supported: Rayleigh proportional damping (C = α M + β K)."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"RAYLEIGH"})
                    .doc("Damping model type. Only RAYLEIGH is supported.")
        );

        // Capture TYPE and then parse one line: alpha, beta
        auto type = std::make_shared<std::string>();
        command.on_enter([type](const fem::dsl::Keys& keys) {
            *type = keys.raw("TYPE");
        });

        command.variant(
            fem::dsl::Variant::make()
                .doc("One data line: α, β coefficients for Rayleigh damping.")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>()
                                .name("RAYLEIGH")
                                .desc("α, β (mass- and stiffness-proportional coefficients).")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([&parser, type](const std::array<fem::Precision, 2>& ab) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "DAMPING must appear inside *LOADCASE.");

                            // Only Transient currently supports DAMPING.
                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                logging::error(*type == "RAYLEIGH", "DAMPING TYPE must be RAYLEIGH.");
                                fem::damping::Rayleigh r{ab[0], ab[1]};
                                lc->set_damping(r);
                                return;
                            }

                            logging::error(false, "DAMPING not supported for loadcase type " + parser.active_loadcase_type());
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
