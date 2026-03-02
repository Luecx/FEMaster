// register_loadcase_inertiarelief.inl — registers INERTIARELIEF within *LOADCASE (LinearStatic only)

#include <stdexcept>

#include "../parser.h"

#include "../../dsl/keyword.h"
#include "../../loadcase/linear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_inertiarelief(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("INERTIARELIEF", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Enable inertia relief for the active linear static load case. "
                    "CONSIDER_POINT_MASSES controls whether POINTMASS features are included.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("CONSIDER_POINT_MASSES")
                    .doc("If true, inertia relief includes all POINTMASS features")
                    .optional("1")
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("INERTIARELIEF must appear inside *LOADCASE");
            }
            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->inertia_relief = true;
                lc->inertia_relief_consider_point_masses = keys.get<bool>("CONSIDER_POINT_MASSES");
                return;
            }
            throw std::runtime_error("INERTIARELIEF is only supported for LINEARSTATIC load cases");
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
