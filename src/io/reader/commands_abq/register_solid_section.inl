/**
 * @file register_solid_section.inl
 * @brief Registers homogeneous Abaqus *SOLID SECTION definitions.
 *
 * Abaqus uses `*SOLID SECTION` for both continuum solids and truss elements.
 * This registration inspects the target FEMaster element set and maps the
 * section either to `model::Model::solid_section()` or, for a pure T3D2-derived
 * truss set, to `model::Model::truss_section()` using the first data value as
 * cross-sectional area.
 *
 * Composite solids, generalized plane-strain reference nodes and section-control
 * options are intentionally unsupported. The optional material orientation is
 * forwarded for continuum solids and ignored for truss sections because the
 * current FEMaster truss section has no material-orientation state.
 *
 * @see model::Model::solid_section
 * @see model::Model::truss_section
 * @see model::T3
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../../model/truss/truss.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers homogeneous Abaqus solid/truss section assignment.
 *
 * The keyword requires `ELSET` and `MATERIAL` and may reference an orientation.
 * The first data value is interpreted according to the target element family:
 * continuum solids ignore the default element attribute, while trusses use it
 * as their cross-sectional area. A target set mixing solids, trusses or other
 * element families is rejected because one FEMaster section cannot represent
 * those different section semantics.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the section definition.
 */
inline void register_solid_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a homogeneous Abaqus solid or truss section to an element set.");

        // Preserve keyword data until the Abaqus section data line determines
        // whether a continuum or truss section must be created
        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET")
                    .required()
                    .doc("Target element set")
                .key("MATERIAL")
                    .required()
                    .doc("Material assigned to the section")
                .key("ORIENTATION")
                    .optional()
                    .doc("Optional continuum-solid material orientation")
        );

        command.on_enter([material, elset, orientation](const fem::io::dsl::Keys& keys) {
            *material    = keys.raw("MATERIAL");
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
        });

        // Abaqus defines one section attribute data line. For ordinary continuum
        // solids the default value is unused here; for T3D2 it is the truss area.
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("ATTRIBUTE").desc("Element-family section attribute")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model, material, elset, orientation](fem::Precision attribute) {
                    if (!model._data->elem_sets.has(*elset)) {
                        throw std::runtime_error("SOLID SECTION element set '" + *elset + "' does not exist");
                    }

                    auto region = model._data->elem_sets.get(*elset);
                    bool has_solid = false;
                    bool has_truss = false;
                    bool has_other = false;

                    // Resolve the section family from the elements already
                    // constructed in the topology pass
                    for (fem::ID id : *region) {
                        auto& element = model._data->elements[id];
                        if (!element) {
                            throw std::runtime_error("SOLID SECTION references an undefined element");
                        }

                        if (element->as<model::T3>() != nullptr) {
                            has_truss = true;
                            continue;
                        }

                        auto structural = element->as<model::StructuralElement>();
                        if (structural != nullptr && structural->is_solid()) {
                            has_solid = true;
                        } else {
                            has_other = true;
                        }
                    }

                    if (region->size() == 0 || has_other || (has_solid && has_truss)) {
                        throw std::runtime_error("SOLID SECTION requires a non-empty pure solid or pure truss element set");
                    }

                    if (has_truss) {
                        model.truss_section(*elset, *material, attribute);
                        return;
                    }

                    model.solid_section(*elset, *material, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
