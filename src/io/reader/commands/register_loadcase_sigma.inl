/**
 * @file register_loadcase_sigma.inl
 * @brief Registers the spectral shift used by linear buckling extraction.
 *
 * The `SIGMA` child command reads a scalar target shift and stores it on the
 * active `LinearBuckling` load case. The value guides the eigenvalue search
 * toward the desired portion of the buckling spectrum without changing the
 * assembled elastic or geometric stiffness operators.
 *
 * Use on any other load-case type is rejected explicitly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../parser.h"

#include "../../../core/logging.h"
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
                    auto* lc = dynamic_cast<loadcase::LinearBuckling*>(parser.active_loadcase());
                    logging::error(lc != nullptr,
                        "SIGMA only valid for LINEARBUCKLING loadcases");
                    lc->sigma = sigma;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
