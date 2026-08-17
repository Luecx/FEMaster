/**
 * @file register_expansion.inl
 * @brief Registers isotropic Abaqus *EXPANSION material data.
 *
 * FEMaster currently stores one scalar thermal-expansion coefficient on a
 * material. The Abaqus isotropic expansion form therefore maps directly onto
 * the existing material state, while orthotropic, transversely isotropic and
 * anisotropic expansion remain unsupported because they require directional
 * coefficients that FEMaster does not currently store.
 *
 * @see material::Material
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <stdexcept>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers constant isotropic Abaqus thermal expansion for the active material.
 *
 * `TYPE=ISO` is the supported Abaqus form. A single constant coefficient is
 * read and forwarded to `Material::set_thermal_expansion()`. Temperature and
 * field-variable dependencies are intentionally unsupported in this initial
 * reader.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model containing the active material.
 */
inline void register_expansion(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("EXPANSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign constant isotropic Abaqus thermal expansion to the active material.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("ISO")
                    .allowed({"ISO", "ISOTROPIC"})
                    .doc("Isotropic thermal-expansion form")
        );

        // Store the single isotropic coefficient directly on the active
        // FEMaster material
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ALPHA").desc("Thermal expansion coefficient")
                )
                .bind([&model](fem::Precision alpha) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("EXPANSION requires an active material context");
                    }
                    material->set_thermal_expansion(alpha);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
