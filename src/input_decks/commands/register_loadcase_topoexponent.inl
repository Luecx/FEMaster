// register_loadcase_topoexponent.inl — registers TOPOEXPONENT for LINEARSTATICTOPO loadcases

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_topoexponent(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPOEXPONENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set penalization exponent for LINEARSTATICTOPO loadcases.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("EXPONENT").desc("Penalization exponent")
                )
                .bind([&parser](fem::Precision exponent) {
                    auto* lc = parser.active_loadcase_as<loadcase::LinearStaticTopo>();
                    if (!lc) {
                        throw std::runtime_error("TOPOEXPONENT only valid for LINEARSTATICTOPO loadcases");
                    }
                    lc->exponent = exponent;
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

