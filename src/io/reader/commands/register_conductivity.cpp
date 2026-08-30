/**
 * @file register_conductivity.cpp
 * @brief Registers isotropic thermal-conductivity material data.
 *
 * The material-scoped `CONDUCTIVITY` command assigns one constant isotropic
 * thermal conductivity to the active material. Thermal elements consume this
 * coefficient when assembling their conductivity matrices and when recovering
 * Fourier heat flux from a solved temperature field.
 *
 * Native syntax:
 *
 * @code
 * *MATERIAL, NAME=STEEL
 * *CONDUCTIVITY
 * 45.0
 * @endcode
 *
 * The initial thermal formulation accepts one strictly positive scalar value;
 * temperature-dependent and anisotropic conductivity are outside the current
 * scope.
 *
 * @see material::Material
 * @see model::ThermalElement
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers constant isotropic conductivity for the active material.
 *
 * Exactly one strictly positive scalar is parsed. The surrounding `MATERIAL`
 * scope selects the material that receives the property through
 * `Material::set_thermal_conductivity()`.
 *
 * @param registry DSL registry receiving the `CONDUCTIVITY` material command.
 * @param model Model owning the currently active material definition.
 */
void register_conductivity(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONDUCTIVITY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign one positive isotropic thermal conductivity to the active material.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("K").desc("Isotropic thermal conductivity")
                )
                .bind([&model](fem::Precision conductivity) {
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr,
                        "CONDUCTIVITY requires an active material context");
                    logging::error(conductivity > Precision(0),
                        "CONDUCTIVITY must be positive");
                    material->set_thermal_conductivity(conductivity);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
