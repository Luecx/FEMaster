/**
 * @file register_pload.inl
 * @brief Registers scalar surface pressure loads.
 *
 * Each `PLOAD` row assigns a scalar pressure to a surface or surface set and may
 * reference a reusable amplitude. The command creates `bc::PLoad` objects in the
 * requested load collector while retaining the geometric target by shared
 * region ownership.
 *
 * Pressure direction follows the finite-element surface normal. Evaluation of
 * normals, surface Jacobians and consistent nodal forces is deferred to load
 * assembly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../reference.h"
#include "../../../bc/load_p.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_pload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        auto amplitude = std::make_shared<std::string>();
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                .key("AMPLITUDE").optional()
        );
        command.on_enter([&model, amplitude](const fem::io::dsl::Keys& keys) {
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<fem::Precision>().name("P")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, amplitude](const std::string& target, fem::Precision value) {
                    model::SurfaceRegion::Ptr region;
                    if (model._data->surface_sets.has(target)) {
                        region = model._data->surface_sets.get(target);
                    } else {
                        region = std::make_shared<model::SurfaceRegion>("INTERNAL");
                        region->add(io::reader::compiled_surface_id(model, target));
                    }

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!amplitude->empty()) {
                        logging::error(model._data->amplitudes.has(*amplitude),
                            "PLOAD: amplitude ", *amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(*amplitude);
                    }

                    auto load = std::make_shared<bc::PLoad>();
                    load->region_    = std::move(region);
                    load->pressure_  = value;
                    load->amplitude_ = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
