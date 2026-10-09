/**
 * @file register_solid_section.cpp
 * @brief Registers solid and truss sections through a shared section command.
 *
 * A section may address solids, trusses or both. Mixed regions are split into
 * registered part-local element subsets so that the concrete section types
 * retain correct element associations during Instance compilation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "../../../model/model.h"
#include "../../../model/truss/truss.h"
#include "../../../section/section_solid.h"
#include "../../../section/section_truss.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

void register_solid_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SOLIDSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART"}));

        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto area        = std::make_shared<Precision>(std::numeric_limits<Precision>::quiet_NaN());

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required()
                .key("ELSET").required()
                .key("ORIENTATION").optional()
        );

        command.on_enter([&model, material, elset, orientation, area](const fem::io::dsl::Keys& keys) {
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "SOLIDSECTION: no active part is available");

            *material    = keys.raw("MATERIAL");
            *elset       = keys.raw("ELSET");
            *orientation = keys.raw("ORIENTATION");
            *area        = std::numeric_limits<Precision>::quiet_NaN();

            logging::error(part->elem_sets.has(*elset),
                "SOLIDSECTION: element set ", *elset, " is not defined in part ", part->name);
            logging::error(model._data->materials.has(*material),
                "SOLIDSECTION: material ", *material, " is not defined");
            logging::error(orientation->empty() || model._data->coordinate_systems.has(*orientation),
                "SOLIDSECTION: coordinate system ", *orientation, " is not defined");
        });

        command.on_exit([&model, material, elset, orientation, area](const fem::io::dsl::Keys&) {
            const auto part   = model._data->parts.get();
            const auto region = part->elem_sets.get(*elset);

            std::vector<ID> solid_ids;
            std::vector<ID> truss_ids;

            for (const ID id : *region) {
                const auto it = part->elements.find(id);
                logging::error(it != part->elements.end() && it->second != nullptr,
                    "SOLIDSECTION: element ", id, " is not defined in part ", part->name);

                if (it->second->as<model::T3>()) {
                    truss_ids.push_back(id);
                    continue;
                }

                const auto* structural = it->second->as<model::StructuralElement>();
                logging::error(structural != nullptr && structural->is_solid(),
                    "SOLIDSECTION: element ", id, " is neither a solid nor a truss");
                solid_ids.push_back(id);
            }

            logging::error(!solid_ids.empty() || !truss_ids.empty(),
                "SOLIDSECTION: element set ", *elset, " is empty");
            logging::error(truss_ids.empty() || (std::isfinite(*area) && *area > Precision(0)),
                "SOLIDSECTION: truss elements require a positive cross-sectional area");
            logging::error(solid_ids.empty() || !truss_ids.empty() || std::isnan(*area),
                "SOLIDSECTION: solid elements do not accept an area data line");

            model::ElementRegion::Ptr solid_region = region;
            model::ElementRegion::Ptr truss_region = region;

            // Model::add_section requires registered regions. Split only mixed
            // element sets; pure sets retain the original user-defined region.
            if (!solid_ids.empty() && !truss_ids.empty()) {
                const std::string prefix = "__SOLIDSECTION_" + std::to_string(part->sections.size());
                const std::string solid_name = prefix + "_SOLID";
                const std::string truss_name = prefix + "_TRUSS";

                logging::error(!part->elem_sets.has(solid_name) && !part->elem_sets.has(truss_name),
                    "SOLIDSECTION: reserved subset names already exist in part ", part->name);

                solid_region = part->elem_sets.activate(solid_name);
                truss_region = part->elem_sets.activate(truss_name);
                for (const ID id : solid_ids) solid_region->add(id);
                for (const ID id : truss_ids) truss_region->add(id);
            }

            const auto assigned_material = model._data->materials.get(*material);

            if (!solid_ids.empty()) {
                auto section = std::make_shared<SolidSection>();
                section->material_    = assigned_material;
                section->region_      = solid_region;
                section->orientation_ = orientation->empty()
                    ? nullptr : model._data->coordinate_systems.get(*orientation);
                model.add_section(std::move(section));
            }

            if (!truss_ids.empty()) {
                model.add_section(std::make_shared<TrussSection>(
                    assigned_material, truss_region, *area));
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("AREA")
                )
                .bind([area](Precision value) {
                    *area = value;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
