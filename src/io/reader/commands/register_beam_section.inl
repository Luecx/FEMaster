/**
 * @file register_beam_section.inl
 * @brief Registers the BEAMSECTION input command.
 *
 * `BEAMSECTION` assigns an existing material and beam profile to an element set
 * in the active semantic Part. Its single data line supplies the non-zero `n1`
 * direction used by beam elements to construct their transverse local basis.
 *
 * `Model::add_section()` retains the resulting `BeamSection` on the active Part.
 * During model compilation the assignment is copied for every Instance, its
 * element region is remapped to assembly identifiers and its spatial direction
 * is rotated with the Instance placement.
 *
 * @see BeamSection
 * @see model::Model::add_section
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../model/model.h"
#include "../../../section/section_beam.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers material, profile and orientation assignments for beam elements.
 *
 * Keyword values and the `n1` data vector are retained across the command's
 * entry, data and exit callbacks. Once the complete command has been consumed,
 * the referenced definitions are resolved against the active Part and model,
 * and the resulting `BeamSection` is registered with the model.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Model providing the active Part and referenced definitions.
 */
inline void register_beam_section(dsl::Registry& registry, model::Model& model) {
    registry.command("BEAMSECTION", [&](dsl::Command& command) {
        // Admit section assignments in the implicit root Part and explicit Parts
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART"}));

        // Expose the command purpose through the registry documentation
        command.doc(
            "Assign a material and beam profile to an element set using the "
            "n1 direction supplied on the following data line."
        );

        // Allocate state shared by the registered entry, data and exit callbacks
        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto profile     = std::make_shared<std::string>();
        auto orientation = std::make_shared<Vec3>(Vec3::Zero());

        // Define the named material, target element set and beam profile
        command.keyword(
            dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required().doc("Material assigned to the beam elements")
                .key("ELSET").required().doc("Element set receiving the section assignment")
                .key("PROFILE").required().doc("Beam profile providing the cross-section properties")
        );

        // Initialize the shared callback state for the current command occurrence
        command.on_enter([material, elset, profile, orientation](const dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
            *profile  = keys.raw("PROFILE");
            orientation->setZero();
        });

        // Resolve all references and create the section after its data line was consumed
        command.on_exit([&model, material, elset, profile, orientation](const dsl::Keys&) {
            // Obtain the semantic Part selected by the surrounding parser scope
            const auto part = model._data->parts.get();

            // Validate the target region, shared definitions and orientation vector
            logging::error(part != nullptr,
                "BEAMSECTION: no active part is available");
            logging::error(part->elem_sets.has(*elset),
                "BEAMSECTION: element set ", *elset, " is not defined in part ", part->name);
            logging::error(model._data->materials.has(*material),
                "BEAMSECTION: material ", *material, " is not defined");
            logging::error(model._data->profiles.has(*profile),
                "BEAMSECTION: profile ", *profile, " is not defined");
            logging::error(orientation->norm() > Precision(0),
                "BEAMSECTION: orientation vector must be non-zero");

            // Assemble the resolved Part-local beam section assignment
            auto section = std::make_shared<BeamSection>();
            section->material_  = model._data->materials.get(*material);
            section->region_    = part->elem_sets.get(*elset);
            section->profile_   = model._data->profiles.get(*profile);
            section->direction_ = *orientation;
            model.add_section(std::move(section));
        });

        // Read the required n1 direction in the active Part coordinate system
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make()
                    .fixed<Precision, 3>().name("N1").desc("Beam-section orientation direction")
                )
                .bind([orientation](const std::array<Precision, 3>& values) {
                    // Retain the parsed direction until the exit callback creates the section
                    *orientation = Vec3{values[0], values[1], values[2]};
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
