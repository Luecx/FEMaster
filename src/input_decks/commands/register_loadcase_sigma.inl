// register_loadcase_sigma.inl — registers SIGMA for LINEARBUCKLING loadcases

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"

namespace fem::input_decks::commands {

inline void register_loadcase_sigma(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("SIGMA", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set eigen-shift parameter for LINEARBUCKLING loadcases.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("SIGMA").desc("Shift parameter")
                )
                .bind([&parser](fem::Precision sigma) {
                    auto* lc = parser.active_loadcase_as<loadcase::LinearBuckling>();
                    if (!lc) {
                        throw std::runtime_error("SIGMA only valid for LINEARBUCKLING loadcases");
                    }
                    lc->sigma = sigma;
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

