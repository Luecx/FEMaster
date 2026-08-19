/**
 * @file register_density.inl
 * @brief Registers material density definitions for FEMaster input readers.
 *
 * The registration attaches a scalar mass density to the material activated by
 * the surrounding `MATERIAL` command. Material storage and active-material
 * state remain owned by `model::Model` and its material repository.
 *
 * The value is validated and retained as constitutive material data for element
 * mass assembly, body loads, inertia relief and dynamic analyses. This command
 * does not construct any mass matrix itself.
 *
 * @see model::Model
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#pragma once

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/registry.h"

namespace fem::io::reader::commands {

/**
 * Registers the `DENSITY` command and assigns its scalar value to the active
 * material.
 *
 * The command is valid only below a material definition. Execution resolves
 * the currently active material, validates that the parser has established the
 * required material context and stores the supplied density on that material.
 *
 * @param registry DSL registry receiving the density command.
 * @param model FEMaster model containing the active material.
 */
inline void register_density(fem::io::dsl::Registry& registry, model::Model& model) {
    // Restrict density definitions to an active material block
    registry.command("DENSITY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign density to the active material.");

        // Parse one scalar density value and apply it to the active material
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("RHO").desc("Density value")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision rho) {
                    // Resolve the material selected by the surrounding command
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr,
                        "DENSITY requires an active material context");
                    material->set_density(rho);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
