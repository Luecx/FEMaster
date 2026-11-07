/**
 * @file register_conductivity.cpp
 * @brief Reads constant isotropic conductivity for uncoupled heat transfer.
 *
 * Abaqus CONDUCTIVITY supports more general material tables. This implementation
 * deliberately accepts only one positive, temperature-independent scalar k.
 * Solid thermal elements use that same k in Fourier's law and in their
 * reference-configuration conduction operator.
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <cmath>

namespace fem::io::reader::commands {

/**
 * Registers a single isotropic heat-conductivity coefficient on MATERIAL.
 *
 * Requiring precisely one data row prevents a tabulated, temperature-dependent
 * conductivity from silently being interpreted as a constant property.
 */
void register_conductivity(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONDUCTIVITY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign constant isotropic heat conductivity to the active material.");

        command.keyword(fem::io::dsl::KeywordSpec::make()
            .key("TYPE").optional("ISO").allowed({"ISO", "ISOTROPIC"})
        );

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("K").desc("Constant isotropic conductivity")
                )
                .bind([&model](Precision conductivity) {
                    logging::error(std::isfinite(conductivity) && conductivity > Precision(0),
                        "CONDUCTIVITY requires a finite, positive coefficient");
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr,
                        "CONDUCTIVITY requires an active material context");
                    material->set_thermal_conductivity(conductivity);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
