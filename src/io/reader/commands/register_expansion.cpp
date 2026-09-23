/**
 * @file register_expansion.cpp
 * @brief Registers isotropic expansion and its legacy keyword alias.
 *
 * The material-scoped `EXPANSION` command and legacy `THERMALEXPANSION` alias
 * use the same constant isotropic grammar in the native and Abaqus readers.
 * Both assign one coefficient to the active `MATERIAL` definition.
 *
 * The stored coefficient is later consumed together with temperature fields by
 * thermal load and constitutive calculations; this file only defines the input
 * grammar and transfers the scalar material property.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <initializer_list>
#include <string>

namespace fem::io::reader::commands {

/**
 * Registers the current and legacy names for one isotropic expansion property.
 *
 * `TYPE=ISO` and `TYPE=ISOTROPIC` are equivalent; omitting TYPE selects ISO.
 * Each command reads one constant alpha with units of inverse temperature.
 * The coefficient is stored on the material active during the definition pass;
 * temperature fields and reference temperatures are applied later by thermal
 * load assembly. Both input dialects use the same grammar and material state.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model containing the active material.
 */
void register_expansion(fem::io::dsl::Registry& registry, model::Model& model) {
    // Register both spellings with identical parameters and scalar data layout.
    for (const char* name : {"EXPANSION", "THERMALEXPANSION"}) {
        registry.command(name, [&, name](fem::io::dsl::Command& command) {
            command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
            command.doc(name == std::string("EXPANSION")
                ? "Assign constant isotropic thermal expansion to the active material."
                : "Legacy alias for EXPANSION; assigns constant isotropic thermal expansion.");

            command.keyword(
                fem::io::dsl::KeywordSpec::make()
                    .key("TYPE")
                        .optional("ISO")
                        .allowed({"ISO", "ISOTROPIC"})
                        .doc("Isotropic thermal-expansion form")
            );

            // Transfer alpha to the active material; the keyword has no load state.
            command.variant(fem::io::dsl::Variant::make()
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1).max(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .one<fem::Precision>().name("ALPHA").desc("Thermal expansion coefficient")
                    )
                    .bind([&model, name](fem::Precision alpha) {
                        auto material = model._data->materials.get();
                        logging::error(material != nullptr,
                            name, " requires an active material context");
                        material->set_thermal_expansion(alpha);
                    })
                )
            );
        });
    }
}

} // namespace fem::io::reader::commands
