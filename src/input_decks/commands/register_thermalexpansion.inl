// register_thermalexpansion.inl — DSL registration for *THERMALEXPANSION

#include <stdexcept>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"

namespace fem::input_decks::commands {

inline void register_thermal_expansion(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("THERMALEXPANSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign thermal expansion coefficient to the active material.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ALPHA").desc("Thermal expansion coefficient")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision alpha) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("THERMALEXPANSION requires an active material context");
                    }
                    material->set_thermal_expansion(alpha);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
