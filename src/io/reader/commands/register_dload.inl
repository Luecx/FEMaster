/**
 * @file register_dload.inl
 * @brief Registers distributed surface tractions.
 *
 * Each `DLOAD` row associates a surface or surface set with a three-component
 * distributed traction. Optional coordinate-system and amplitude references are
 * retained by `bc::DLoad`, and the completed object is inserted into the named
 * load collector.
 *
 * Surface integration, interpolation to element nodes and physical-area scaling
 * remain responsibilities of the load and surface implementations.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../bc/load_d.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_dload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
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
                    .fixed<fem::Precision, 3>().name("LOAD")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, orientation, amplitude](const std::string& target,
                                                       const std::array<fem::Precision, 3>& values) {
                    model::SurfaceRegion::Ptr region;
                    if (model._data->surface_sets.has(target)) {
                        region = model._data->surface_sets.get(target);
                    } else {
                        region = std::make_shared<model::SurfaceRegion>("INTERNAL");
                        region->add(io::reader::compiled_surface_id(model, target));
                    }

                    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
                    if (!orientation->empty()) {
                        logging::error(model._data->coordinate_systems.has(*orientation),
                            "DLOAD: coordinate system ", *orientation, " does not exist");
                        orientation_ptr = model._data->coordinate_systems.get(*orientation);
                    }

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!amplitude->empty()) {
                        logging::error(model._data->amplitudes.has(*amplitude),
                            "DLOAD: amplitude ", *amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(*amplitude);
                    }

                    auto load = std::make_shared<bc::DLoad>();
                    load->region_      = std::move(region);
                    load->values_      = Vec3{values[0], values[1], values[2]};
                    load->orientation_ = std::move(orientation_ptr);
                    load->amplitude_   = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
