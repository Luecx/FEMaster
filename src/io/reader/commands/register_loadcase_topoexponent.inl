/**
 * @file register_loadcase_topoexponent.inl
 * @brief Registers the density penalization exponent for topology statics.
 *
 * The `TOPOEXPONENT` child command reads the scalar exponent used by an active
 * `LinearStaticTopo` analysis. It configures the SIMP-like material interpolation
 * applied to element densities when stiffness and recovered quantities are
 * evaluated.
 *
 * The command validates the active formulation and stores the parameter; the
 * actual penalization remains part of the topology-aware load case.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_static_topo.h"

namespace fem::io::reader::commands {

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
                    auto* lc = dynamic_cast<loadcase::LinearStaticTopo*>(parser.active_loadcase());
                    logging::error(lc != nullptr,
                        "TOPOEXPONENT only valid for LINEARSTATICTOPO loadcases");
                    lc->exponent = exponent;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
