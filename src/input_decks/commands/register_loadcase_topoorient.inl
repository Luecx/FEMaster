// register_loadcase_topoorient.inl — registers TOPOORIENT for LINEARSTATICTOPO loadcases

#include <array>
#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_topoorient(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPOORIENT", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set orientation angles for LINEARSTATICTOPO loadcases.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("EID").desc("Element id")
                    .fixed<fem::Precision, 3>().name("ANGLE").desc("Orientation angles")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&parser](fem::ID id, const std::array<fem::Precision, 3>& angles) {
                    auto* lc = parser.active_loadcase_as<loadcase::LinearStaticTopo>();
                    if (!lc) {
                        throw std::runtime_error("TOPOORIENT only valid for LINEARSTATICTOPO loadcases");
                    }
                    lc->orientation(id, 0) = angles[0];
                    lc->orientation(id, 1) = angles[1];
                    lc->orientation(id, 2) = angles[2];
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

