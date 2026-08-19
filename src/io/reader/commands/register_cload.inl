/**
 * @file register_cload.inl
 * @brief Registers concentrated nodal loads.
 *
 * Each `CLOAD` row targets a node or node set and supplies six generalized force
 * and moment components. Optional coordinate-system and amplitude references
 * are resolved once and stored with the resulting `bc::CLoad` in the selected
 * load collector.
 *
 * Instance-qualified node references are translated through compiled assembly
 * maps. Coordinate transformation and time-dependent scaling are applied by the
 * load during analysis rather than baked into the parsed values.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../bc/load_c.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_cload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        auto orientation = std::make_shared<std::string>();
        auto amplitude   = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                .key("ORIENTATION").optional()
                .key("AMPLITUDE").optional()
        );

        command.on_enter([&model, orientation, amplitude](const fem::io::dsl::Keys& keys) {
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *amplitude   = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .fixed<fem::Precision, 6>().name("LOAD")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, orientation, amplitude](const std::string& target,
                                                       const std::array<fem::Precision, 6>& values) {
                    model::NodeRegion::Ptr region;
                    if (model._data->node_sets.has(target)) {
                        region = model._data->node_sets.get(target);
                    } else {
                        region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(io::reader::compiled_node_id(model, target));
                    }

                    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
                    if (!orientation->empty()) {
                        logging::error(model._data->coordinate_systems.has(*orientation),
                            "CLOAD: coordinate system ", *orientation, " does not exist");
                        orientation_ptr = model._data->coordinate_systems.get(*orientation);
                    }

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!amplitude->empty()) {
                        logging::error(model._data->amplitudes.has(*amplitude),
                            "CLOAD: amplitude ", *amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(*amplitude);
                    }

                    auto load = std::make_shared<bc::CLoad>();
                    load->region_      = std::move(region);
                    load->orientation_ = std::move(orientation_ptr);
                    load->amplitude_   = std::move(amplitude_ptr);
                    load->values_ << values[0], values[1], values[2], values[3], values[4], values[5];
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
