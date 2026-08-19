/**
 * @file register_loadcase_supports.inl
 * @brief Registers support and constraint collectors for active load cases.
 *
 * The `SUPPORTS` child command appends non-empty collector names to every
 * supported structural analysis, including static, buckling, eigenfrequency,
 * transient, harmonic and nonlinear formulations. Input order is retained so
 * downstream constraint construction remains deterministic.
 *
 * Resolution of support, tie and coupling collectors and construction of the
 * resulting constraint equations remain responsibilities of the load case and
 * model constraint subsystem.
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
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_static_topo.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_supports(fem::io::dsl::Registry& registry, Parser& parser) {
    const auto append_tokens = [](const std::array<std::string, 16>& tokens, std::vector<std::string>& out) {
        for (const auto& token : tokens) {
            if (!token.empty()) out.push_back(token);
        }
    };

    registry.command("SUPPORTS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Assign support collectors to the active loadcase.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<std::string, 16>().name("SUPP").desc("Support collector names")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser, append_tokens](const std::array<std::string, 16>& names) {
                    auto* base = parser.active_loadcase();
                    logging::error(base != nullptr,
                        "SUPPORTS must appear inside *LOADCASE");

                    if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearStaticTopo*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }

                    logging::error(false,
                        "SUPPORTS not supported for loadcase type ", base->type_name());
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
