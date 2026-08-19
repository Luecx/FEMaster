/**
 * @file register_loadcase_inertiarelief.inl
 * @brief Registers inertia-relief settings for linear static load cases.
 *
 * `INERTIARELIEF` enables rigid-body force and moment balancing on the active
 * `LinearStatic` analysis. The optional `CONSIDER_POINT_MASSES` setting controls
 * whether concentrated mass features participate in the mass, inertia and
 * compensating-load calculation.
 *
 * The command only configures the analysis. Construction of the balancing
 * inertial load remains the responsibility of the linear static solver.
 * Derived linear-static formulations inherit the same setting.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../parser.h"

#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../loadcase/linear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_inertiarelief(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("INERTIARELIEF", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable inertia relief for the active linear static load case. "
                    "CONSIDER_POINT_MASSES controls whether POINTMASS features are included.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("CONSIDER_POINT_MASSES")
                    .doc("If true, inertia relief includes all POINTMASS features")
                    .optional("1")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            logging::error(base != nullptr,
                "INERTIARELIEF must appear inside *LOADCASE");

            auto* lc = dynamic_cast<loadcase::LinearStatic*>(base);
            logging::error(lc != nullptr,
                "INERTIARELIEF is only supported for LINEARSTATIC load cases");
            lc->inertia_relief = true;
            lc->inertia_relief_consider_point_masses = keys.get<bool>("CONSIDER_POINT_MASSES");
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
