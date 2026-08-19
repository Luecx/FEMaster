/**
 * @file register_loadcase_damping.inl
 * @brief Registers Rayleigh damping for transient and harmonic analyses.
 *
 * The `DAMPING` child command reads the mass- and stiffness-proportional
 * coefficients of `C = alpha M + beta K`. It constructs the common
 * `RayleighDamping` representation and assigns it to an active linear transient
 * or direct harmonic load case.
 *
 * Assembly of the physical or reduced damping operator remains part of the
 * concrete analysis because its configuration and constraint transformation
 * depend on the selected solution procedure.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../core/types_eig.h"

#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/tools/rayleigh_damping.h"

namespace fem::io::reader::commands {

inline void register_loadcase_damping(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("DAMPING", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Assign Rayleigh proportional damping C = alpha M + beta K.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"RAYLEIGH"})
                    .doc("Damping model type. Only RAYLEIGH is supported.")
        );

        auto type = std::make_shared<std::string>();
        command.on_enter([type](const fem::io::dsl::Keys& keys) {
            *type = keys.raw("TYPE");
        });

        command.variant(
            fem::io::dsl::Variant::make()
                .doc("One data line: alpha, beta coefficients for Rayleigh damping.")
                .segment(
                    fem::io::dsl::Segment::make()
                        .range(fem::io::dsl::LineRange{}.min(1).max(1))
                        .pattern(
                            fem::io::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>()
                                .name("RAYLEIGH")
                                .desc("alpha, beta mass- and stiffness-proportional coefficients.")
                                .on_missing(fem::Precision{0})
                                .on_empty(fem::Precision{0})
                        )
                        .bind([&parser, type](const std::array<fem::Precision, 2>& ab) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "DAMPING must appear inside *LOADCASE.");
                            logging::error(*type == "RAYLEIGH", "DAMPING TYPE must be RAYLEIGH.");

                            const fem::loadcase::tools::RayleighDamping damping{ab[0], ab[1]};

                            if (auto* lc = dynamic_cast<fem::loadcase::Transient*>(base)) {
                                lc->set_damping(damping);
                                return;
                            }

                            if (auto* lc = dynamic_cast<fem::loadcase::LinearHarmonic*>(base)) {
                                lc->set_damping(damping);
                                return;
                            }

                            logging::error(false,
                                           "DAMPING not supported for loadcase type " +
                                           base->type_name());
                        })
                )
        );
    });
}

} // namespace fem::io::reader::commands
