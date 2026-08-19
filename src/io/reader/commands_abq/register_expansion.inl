/**
 * @file register_expansion.inl
 * @brief Registers isotropic Abaqus *EXPANSION material data.
 *
 * Constant isotropic thermal expansion is mapped to the scalar thermal-expansion
 * coefficient stored by FEMaster materials. Direction-dependent and
 * temperature-/field-dependent expansion forms are outside the supported syntax.
 *
 * The command executes while the surrounding Abaqus material is active during
 * the definition pass. Thermal strains and equivalent forces are evaluated
 * later from this coefficient, reference temperatures and nodal temperature
 * fields; no analysis state is constructed here.
 *
 * @see material::Material
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers constant isotropic Abaqus thermal expansion for the active material.
 *
 * `TYPE=ISO` and `TYPE=ISOTROPIC` read one constant coefficient and assign it to
 * the material selected by the enclosing `*MATERIAL` block.
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

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ALPHA").desc("Thermal expansion coefficient")
                )
                .bind([&model](fem::Precision alpha) {
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr,
                        "EXPANSION requires an active material context");
                    material->set_thermal_expansion(alpha);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
