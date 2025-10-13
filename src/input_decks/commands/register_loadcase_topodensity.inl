// register_loadcase_topodensity.inl — registers TOPODENSITY for LINEARSTATICTOPO loadcases

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_topodensity(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPODENSITY", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set element densities for LINEARSTATICTOPO loadcases.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("EID").desc("Element id")
                    .one<fem::Precision>().name("DENS").desc("Element density")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&parser](fem::ID id, fem::Precision density) {
                    auto* lc = parser.active_loadcase_as<loadcase::LinearStaticTopo>();
                    if (!lc) {
                        throw std::runtime_error("TOPODENSITY only valid for LINEARSTATICTOPO loadcases");
                    }
                    lc->density(id) = density;
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
