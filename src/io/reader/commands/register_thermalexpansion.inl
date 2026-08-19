/**
 * @file register_thermalexpansion.inl
 * @brief Registers isotropic thermal-expansion material data.
 *
 * The material-scoped `THERMALEXPANSION` command reads one coefficient of
 * thermal expansion and assigns it to the material active in the surrounding
 * `MATERIAL` definition. Validation ensures that the command cannot silently
 * modify an absent material context.
 *
 * The stored coefficient is later consumed together with temperature fields by
 * thermal load and constitutive calculations; this file only defines the input
 * grammar and transfers the scalar material property.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../dsl/condition.h"

namespace fem::io::reader::commands {

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
                    logging::error(material != nullptr,
                        "THERMALEXPANSION requires an active material context");
                    material->set_thermal_expansion(alpha);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
