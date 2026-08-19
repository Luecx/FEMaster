/**
 * @file register_shell_section.inl
 * @brief Registers integrated and ABD shell sections.
 *
 * `SHELLSECTION` supports material-integrated thickness definitions and direct
 * membrane-bending `ABD` constitutive data. It resolves the target element set,
 * optional material orientation and coordinate-system axis before creating the
 * corresponding shell-section implementation.
 *
 * Integrated sections defer through-thickness constitutive evaluation to their
 * material, whereas ABD sections store the supplied generalized stiffness and
 * inertia terms directly. Both are assigned before model compilation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>

#include "../../../model/model.h"
#include "../../../section/section_shell_abd.h"
#include "../../../section/section_shell_integrated.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_shell_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto thickness   = std::make_shared<Precision>(Precision(1));
        auto csys_axis   = std::make_shared<Index>(0);

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("INTEGRATED").allowed({"INTEGRATED", "ABD"})
                .key("MATERIAL").alternative("MAT").optional()
                .key("ELSET").required()
                .key("THICKNESS").optional("1.0")
                .key("ORIENTATION").optional()
                .key("CSYSAXIS").optional("1").allowed({"1", "2", "3"})
        );
        command.on_enter([material, elset, orientation, thickness, csys_axis](const fem::io::dsl::Keys& keys) {
            *material    = keys.has("MATERIAL") ? keys.raw("MATERIAL") : std::string{};
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *thickness   = keys.get<Precision>("THICKNESS");
            *csys_axis   = static_cast<Index>(keys.get<int>("CSYSAXIS") - 1);
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"INTEGRATED"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("THICKNESS")
                        .on_missing(Precision{1}).on_empty(Precision{1})
                )
                .bind([&model, material, elset, orientation, csys_axis](Precision local_thickness) {
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SHELLSECTION: no active part is available");
                    logging::error(part->elem_sets.has(*elset),
                        "SHELLSECTION: element set ", *elset, " is not defined in part ", part->name);
                    logging::error(!material->empty() && model._data->materials.has(*material),
                        "SHELLSECTION: integrated section requires a defined material");
                    logging::error(orientation->empty() || model._data->coordinate_systems.has(*orientation),
                        "SHELLSECTION: coordinate system ", *orientation, " is not defined");

                    model.add_section(std::make_shared<IntegratedShellSection>(
                        model._data->materials.get(*material),
                        part->elem_sets.get(*elset),
                        local_thickness,
                        orientation->empty() ? nullptr : model._data->coordinate_systems.get(*orientation),
                        *csys_axis
                    ));
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ABD"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(8))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<Precision, 40>().name("DATA")
                )
                .bind([&model, material, elset, orientation, thickness, csys_axis](const std::array<Precision, 40>& values) {
                    const auto part = model._data->parts.get();
                    logging::error(part != nullptr,
                        "SHELLSECTION: no active part is available");
                    logging::error(part->elem_sets.has(*elset),
                        "SHELLSECTION: element set ", *elset, " is not defined in part ", part->name);
                    logging::error(material->empty() || model._data->materials.has(*material),
                        "SHELLSECTION: material ", *material, " is not defined");
                    logging::error(orientation->empty() || model._data->coordinate_systems.has(*orientation),
                        "SHELLSECTION: coordinate system ", *orientation, " is not defined");

                    StaticMatrix<6, 6> abd;
                    StaticMatrix<2, 2> shear;
                    for (Index i = 0; i < 6; ++i) {
                        for (Index j = 0; j < 6; ++j) abd(i, j) = values[6 * i + j];
                    }
                    for (Index i = 0; i < 2; ++i) {
                        for (Index j = 0; j < 2; ++j) shear(i, j) = values[36 + 2 * i + j];
                    }

                    model.add_section(std::make_shared<ABDShellSection>(
                        material->empty() ? nullptr : model._data->materials.get(*material),
                        part->elem_sets.get(*elset),
                        *thickness,
                        abd,
                        shear,
                        orientation->empty() ? nullptr : model._data->coordinate_systems.get(*orientation),
                        *csys_axis
                    ));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
