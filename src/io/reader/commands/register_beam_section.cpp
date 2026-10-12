/**
 * @file register_beam_section.cpp
 * @brief Registers native and Abaqus general beam section commands.
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

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include "../../../material/material.h"
#include "../../../material/isotropic_elasticity.h"
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
void register_beam_section(dsl::Registry& registry, model::Model& model) {
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

    // Abaqus generalized beam section with directly specified geometric properties
    registry.command("BEAM GENERAL SECTION", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "PART"}));
        command.doc("Define A, I11, I12, I22 and J, using MATERIAL or direct E, G and alpha.");

        auto elset       = std::make_shared<std::string>();
        auto material    = std::make_shared<std::string>();
        auto values      = std::make_shared<std::array<Precision, 5>>();
        auto direction   = std::make_shared<Vec3>(Vec3{0, 0, -1});
        auto elastic     = std::make_shared<std::array<Precision, 3>>();
        auto density     = std::make_shared<Precision>(Precision(-1));

        command.keyword(dsl::KeywordSpec::make()
            .key("ELSET").required()
            .key("SECTION").optional("GENERAL").allowed({"GENERAL"})
            .key("MATERIAL").optional()
            .key("DENSITY").optional());

        // Keep all section-local data separate between command occurrences
        command.on_enter([=](const dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
            *material = keys.raw("MATERIAL");
            values->fill(Precision(0));
            elastic->fill(Precision(0));
            *direction = Vec3{0, 0, -1};
            *density = keys.has("DENSITY") ? keys.get<Precision>("DENSITY") : Precision(-1);
            logging::error(*density >= Precision(0) || *density == Precision(-1),
                "BEAM GENERAL SECTION: DENSITY must be non-negative");
        });

        // Resolve the generated profile and its material into the standard BeamSection
        command.on_exit([=, &model](const dsl::Keys&) {
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "BEAM GENERAL SECTION: no active part");
            logging::error(part->elem_sets.has(*elset),
                "BEAM GENERAL SECTION: undefined ELSET ", *elset);
            const auto& v = *values;
            logging::error(v[0] > Precision(0) && v[1] >= Precision(0) &&
                           v[3] >= Precision(0) && v[4] >= Precision(0) &&
                           v[1] * v[3] >= v[2] * v[2],
                "BEAM GENERAL SECTION: invalid A, I11, I12, I22 or J");
            logging::error(direction->norm() > Precision(0),
                "BEAM GENERAL SECTION: n1 must be non-zero");

            material::MaterialPtr assigned;
            if (!material->empty()) {
                logging::error(model._data->materials.has(*material),
                    "BEAM GENERAL SECTION: unknown MATERIAL ", *material);
                assigned = model._data->materials.get(*material);
            } else {
                const Precision E = (*elastic)[0];
                const Precision G = (*elastic)[1];
                logging::error(E > Precision(0) && G > Precision(0),
                    "BEAM GENERAL SECTION: direct properties require positive E and G");
                const Precision nu = E / (Precision(2) * G) - Precision(1);
                logging::error(nu > Precision(-1) && nu < Precision(0.5),
                    "BEAM GENERAL SECTION: E and G imply invalid Poisson ratio");
                assigned = std::make_shared<material::Material>("__BEAM_GENERAL_" + *elset);
                assigned->set_elasticity<material::IsotropicElasticity>(E, nu);
                assigned->set_thermal_expansion((*elastic)[2]);
            }

            // A density override is section-local and must not change the named material
            if (*density >= Precision(0)) {
                if (!material->empty()) {
                    const auto* iso = assigned->elasticity()->as<material::IsotropicElasticity>();
                    logging::error(iso != nullptr,
                        "BEAM GENERAL SECTION: DENSITY override requires isotropic material");
                    auto copy = std::make_shared<material::Material>("__BEAM_GENERAL_" + *elset);
                    copy->set_elasticity<material::IsotropicElasticity>(iso->youngs, iso->poisson);
                    if (assigned->has_thermal_expansion()) {
                        copy->set_thermal_expansion(assigned->get_thermal_expansion());
                        copy->set_thermal_zero_temperature(assigned->get_thermal_zero_temperature());
                    }
                    assigned = std::move(copy);
                }
                assigned->set_density(*density);
            }

            auto section = std::make_shared<BeamSection>();
            section->region_ = part->elem_sets.get(*elset);
            section->material_ = std::move(assigned);
            section->profile_ = std::make_shared<Profile>(
                "__BEAM_GENERAL_" + *elset, v[0], v[1], v[3], v[4], v[2]);
            section->direction_ = *direction;
            model.add_section(std::move(section));
        });

        auto geometry = [&]() {
            return dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make().fixed<Precision, 5>().name("A_I11_I12_I22_J"))
                .bind([values](const std::array<Precision, 5>& data) {
                    *values = data;
                });
        };

        auto orientation = [&]() {
            return dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make().fixed<Precision, 3>().name("N1"))
                .bind([direction](const std::array<Precision, 3>& n1) {
                    *direction = Vec3{n1[0], n1[1], n1[2]};
                });
        };

        // With MATERIAL, Abaqus permits the n1 direction to be omitted
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::key_present("MATERIAL"))
            .segment(geometry())
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(0).max(1))
                .pattern(dsl::Pattern::make().fixed<Precision, 3>().name("N1"))
                .bind([direction](const std::array<Precision, 3>& n1) {
                    *direction = Vec3{n1[0], n1[1], n1[2]};
                })));

        // Without MATERIAL, Abaqus requires a third line for E, G and optional alpha
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::negate(dsl::Condition::key_present("MATERIAL")))
            .segment(geometry())
            .segment(orientation())
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1).max(1))
                .pattern(dsl::Pattern::make()
                    .one<Precision>().name("E")
                    .one<Precision>().name("G")
                    .one<Precision>().name("ALPHA")
                        .on_missing(Precision(0)).on_empty(Precision(0)))
                .bind([elastic](Precision E, Precision G, Precision alpha) {
                    *elastic = {E, G, alpha};
                })));
    });

}

} // namespace fem::io::reader::commands
