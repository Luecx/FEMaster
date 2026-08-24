/**
 * @file register_mass.inl
 * @brief Registers isotropic Abaqus mass properties for MASS elements.
 *
 * Abaqus separates one-node `MASS` topology from its scalar property. FEMaster
 * keeps that separation: `model::PointElement` represents only the element id
 * and connectivity, while `*MASS` resolves the point-element connectivity and
 * creates ordinary nodal `feature::PointMass` objects.
 *
 * Part/root MASS records create part-local PointMass features immediately. Their
 * node regions remain in the owning Part's local namespace until
 * `Model::compile()` copies and remaps them for every Instance. Assembly-level
 * MASS records are applied directly against compiled ELSETs during the assembly
 * pass.
 *
 * Only the standard isotropic Abaqus form is supported here. Anisotropic mass,
 * orientation and damping parameters remain outside this translation.
 *
 * @see model::PointElement
 * @see feature::PointMass
 * @see model::Part
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../parser_abq.h"
#include "../../../feature/point_mass.h"
#include "../../../model/element/point.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers `*MASS, ELSET=...` and translates the scalar property to PointMass.
 *
 * Before compilation, non-assembly records resolve the active Part and create
 * one part-local PointMass per referenced PointElement using only its node
 * connectivity. During compilation those features are copied once per Instance.
 * After compilation only assembly-level records execute and create PointMass
 * features directly in assembly space from the compiled PointElements.
 *
 * @param registry Parser registry receiving the MASS command.
 * @param parser Abaqus parser owning the model and assembly-level helper.
 * @param assembly_scope Shared parser flag identifying assembly-level records.
 */
inline void register_mass(fem::io::dsl::Registry& registry,
                          ParserAbq& parser,
                          std::shared_ptr<bool> assembly_scope) {
    registry.command("MASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.doc("Assign isotropic mass to Abaqus MASS elements in an ELSET.");

        auto elset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Element set containing MASS elements")
                .key("TYPE").optional("ISOTROPIC").allowed({"ISOTROPIC"})
                    .doc("Only isotropic point mass is currently supported")
        );

        command.on_enter([elset](const fem::io::dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("MASS").desc("Isotropic mass magnitude")
                )
                .bind([&parser, assembly_scope, elset](Precision mass) {
                    auto& model = parser.model();

                    // Part/root properties belong to the semantic Part itself.
                    // Assembly records are replayed after compilation instead.
                    if (!model._data->compiled) {
                        if (*assembly_scope) return;

                        const auto part = model._data->parts.get();
                        logging::error(part != nullptr,
                            "MASS: no active part is available");
                        logging::error(part->elem_sets.has(*elset),
                            "MASS: element set ", *elset, " does not exist in part ", part->name);

                        const auto region = part->elem_sets.get(*elset);
                        logging::error(region != nullptr && region->size() > 0,
                            "MASS: element set ", *elset, " is empty in part ", part->name);

                        for (const ID id : *region) {
                            const auto it = part->elements.find(id);
                            logging::error(it != part->elements.end() && it->second != nullptr,
                                "MASS: element ", id, " does not exist in part ", part->name);

                            const auto* point = it->second->as<model::PointElement>();
                            logging::error(point != nullptr,
                                "MASS: element set ", *elset, " contains non-MASS element ", id);

                            auto nodes = std::make_shared<model::NodeRegion>(
                                "__ABQ_MASS_" + std::to_string(id));
                            nodes->add(point->nodes()[0]);

                            auto point_mass = std::make_shared<feature::PointMass>();
                            point_mass->region_ = std::move(nodes);
                            point_mass->mass_   = mass;
                            part->point_masses.push_back(std::move(point_mass));
                        }
                        return;
                    }

                    // Part/root records were already compiled from Part storage.
                    // Only assembly-level MASS records execute in this pass.
                    if (!*assembly_scope) return;

                    logging::error(model._data->elem_sets.has(*elset),
                        "MASS: compiled element set ", *elset, " does not exist");
                    const auto region = model._data->elem_sets.get(*elset);
                    logging::error(region != nullptr && region->size() > 0,
                        "MASS: compiled element set ", *elset, " is empty");

                    for (const ID id : *region) {
                        parser.add_abaqus_mass_feature(id, mass);
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
