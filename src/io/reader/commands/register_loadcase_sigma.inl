// register_loadcase_sigma.inl — registers SIGMA for LINEARBUCKLING loadcases

#include <stdexcept>

#include "../parser.h"

#include "../../../loadcase/linear_buckling.h"

namespace fem::io::reader::commands {

inline void register_loadcase_sigma(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("SIGMA", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set eigen-shift parameter for LINEARBUCKLING loadcases.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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

} // namespace fem::io::reader::commands

