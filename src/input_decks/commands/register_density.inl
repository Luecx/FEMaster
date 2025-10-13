// register_density.inl — DSL registration for *DENSITY

#include <stdexcept>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"

namespace fem::input_decks::commands {

inline void register_density(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("DENSITY", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign density to the active material.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("RHO").desc("Density value")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision rho) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("DENSITY requires an active material context");
                    }
                    material->set_density(rho);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
