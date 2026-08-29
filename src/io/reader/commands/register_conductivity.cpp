/**
 * @file register_conductivity.cpp
 * @brief Registers isotropic thermal-conductivity material data.
 *
 * The material-scoped `CONDUCTIVITY` command reads one positive isotropic
 * conductivity and stores it on the material activated by the surrounding
 * `MATERIAL` definition. Solid thermal elements consume this coefficient when
 * assembling their conductivity matrix and recovering conductive heat flux.
 *
 * The command only defines constitutive input. Global thermal assembly remains
 * the responsibility of `model::Model`, while steady-state solution control
 * belongs to `loadcase::SteadyStateThermal`.
 *
 * @see material::Material
 * @see model::ThermalElement
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"

namespace fem::io::reader::commands {

/**
 * Registers constant isotropic conductivity for the active material.
 *
 * Exactly one strictly positive scalar is accepted. A positive coefficient is
 * required by the initial steady-state thermal formulation because zero-valued
 * elements would create disconnected algebraic temperature regions.
 *
 * @param registry DSL registry receiving the material child command.
 * @param model Model owning the currently active material definition.
 */
void register_conductivity(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONDUCTIVITY", [&](fem::io::dsl::Command& command) {
        // Conductivity belongs to the material selected by the parent scope
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign one positive isotropic thermal conductivity to the active material.");

        // Parse and store the constant conductivity coefficient
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
