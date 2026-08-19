/**
 * @file register_loadcase_loads.inl
 * @brief Registers load-collector selection for active load cases.
 *
 * The `LOADS` child command reads one or more load-collector names and appends
 * every non-empty token to the active analysis. It supports the static,
 * buckling, transient, harmonic and nonlinear load cases that assemble external
 * forces from named model collectors.
 *
 * Collector existence and formulation-specific load assembly remain load-case
 * responsibilities. This registration layer validates the active analysis type
 * and preserves the order in which collector names appear in the deck.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include <array>
#include <string>
#include <vector>

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_loads(fem::io::dsl::Registry& registry, Parser& parser) {
    const auto append_tokens = [](const std::array<std::string, 16>& tokens, std::vector<std::string>& out) {
        for (const auto& token : tokens) {
            if (!token.empty()) out.push_back(token);
        }
    };

    registry.command("LOADS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Assign load collectors to the active loadcase.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 16>().name("LOAD").desc("Load collector names")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser, append_tokens](const std::array<std::string, 16>& names) {
                    auto* base = parser.active_loadcase();
                    logging::error(base != nullptr,
                        "LOADS must appear inside *LOADCASE");

                    if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    logging::error(false,
                        "LOADS not supported for loadcase type ", base->type_name());
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
