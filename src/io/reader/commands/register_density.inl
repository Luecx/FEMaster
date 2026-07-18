// register_density.inl — DSL registration for *DENSITY

#include <stdexcept>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"

namespace fem::io::reader::commands {

inline void register_density(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("DENSITY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign density to the active material.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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

} // namespace fem::io::reader::commands
