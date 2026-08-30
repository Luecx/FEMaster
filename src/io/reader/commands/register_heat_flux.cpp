/**
 * @file register_heat_flux.cpp
 * @brief Registers prescribed surface heat fluxes.
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../bc/neumann/heat_flux.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

void register_heat_flux(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("HEATFLUX", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Prescribe a scalar surface heat flux; positive values add heat to the model.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->thermal_load_cols.activate(keys.raw("LOAD_COLLECTOR"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled surface or surface-set reference")
                    .one<fem::Precision>().name("FLUX").desc("Heat flux entering the model")
                )
                .bind([&model](const std::string& target, fem::Precision flux) {
                    model::SurfaceRegion::Ptr region;
                    if (model._data->surface_sets.has(target)) {
                        region = model._data->surface_sets.get(target);
                    } else {
                        region = std::make_shared<model::SurfaceRegion>("INTERNAL");
                        region->add(model.compiled_surface_id(target));
                    }

                    auto load = std::make_shared<bc::HeatFlux>();
                    load->region_ = std::move(region);
                    load->flux_   = flux;

                    logging::error(model._data->thermal_load_cols.get() != nullptr,
                        "HEATFLUX: no thermal load collector is active");
                    model._data->thermal_load_cols.get()->add(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
