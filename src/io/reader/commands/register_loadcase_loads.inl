// register_loadcase_loads.inl — registers LOADS within *LOADCASE

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "../parser.h"

#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_static.h"
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
                .bind([&parser, &append_tokens](const std::array<std::string, 16>& names) {
                    auto* base = parser.active_loadcase();
                    if (!base) {
                        throw std::runtime_error("LOADS must appear inside *LOADCASE");
                    }

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
                    if (auto* lc = dynamic_cast<loadcase::Transient*>(base)) {
                        append_tokens(names, lc->loads);
                        return;
                    }
                    throw std::runtime_error("LOADS not supported for loadcase type " + parser.active_loadcase_type());
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
