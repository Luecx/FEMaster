// register_loadcase_supports.inl — registers SUPPORTS within *LOADCASE

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_eigenfreq.h"
#include "../../loadcase/linear_static.h"
#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_supports(fem::dsl::Registry& registry, Parser& parser) {
    const auto append_tokens = [](const std::array<std::string, 16>& tokens, std::vector<std::string>& out) {
        for (const auto& token : tokens) {
            if (!token.empty()) out.push_back(token);
        }
    };

    registry.command("SUPPORTS", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Assign support collectors to the active loadcase.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<std::string, 16>().name("SUPP").desc("Support collector names")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser, &append_tokens](const std::array<std::string, 16>& names) {
                    auto* base = parser.active_loadcase();
                    if (!base) {
                        throw std::runtime_error("SUPPORTS must appear inside *LOADCASE");
                    }

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
                    if (auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                        append_tokens(names, lc->supps);
                        return;
                    }

                    throw std::runtime_error("SUPPORTS not supported for loadcase type " + parser.active_loadcase_type());
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
