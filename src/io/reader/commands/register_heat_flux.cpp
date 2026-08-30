/**
 * @file register_heat_flux.cpp
 * @brief Registers prescribed thermal surface heat flux.
 *
 * `HEATFLUX` defines a scalar heat flux per unit surface area on a compiled
 * surface or surface set. Positive values denote heat entering the model. The
 * condition is stored in the named common load collector as `bc::HeatFlux` and
 * is later integrated consistently to nodal heat-flow contributions by the
 * thermal analysis.
 *
 * Native syntax:
 *
 * @code
 * *HEATFLUX, LOAD_COLLECTOR=THERMAL_LOADS
 * SURFACE_HOT, 2500.0
 * @endcode
 *
 * An optional `AMPLITUDE` scales the prescribed flux at analysis time. The
 * scalar value from each data row is written directly to `HeatFlux::heat_flux_`;
 * no vector direction or coordinate system is involved.
 *
 * @see bc::HeatFlux
 * @see bc::Neumann
 * @see bc::LoadCollector
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "../../../bc/neumann/heat_flux.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers uniform surface heat-flux loads in the native DSL.
 *
 * Entering the command activates one common load collector and resolves the
 * optional amplitude once for all data rows. Each row resolves `TARGET` first
 * as a compiled surface set and otherwise as one compiled surface reference.
 * `VALUE` is the signed prescribed heat flux and is assigned directly to the
 * created `bc::HeatFlux::heat_flux_` member.
 *
 * @param registry DSL registry receiving the `HEATFLUX` command.
 * @param model Compiled model providing surfaces, amplitudes and load storage.
 */
void register_heat_flux(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("HEATFLUX", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Apply signed heat flux per unit area to compiled surfaces; positive values enter the model.");

        auto amplitude = std::make_shared<std::string>();
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .alternative("LOAD COLLECTOR")
                    .alternative("NAME")
                    .required()
                    .doc("Load collector receiving the heat-flux conditions")
                .key("AMPLITUDE").optional().doc("Optional scalar amplitude")
        );
        command.on_enter([&model, amplitude](const fem::io::dsl::Keys& keys) {
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled surface or surface-set reference")
                    .one<fem::Precision>().name("VALUE").desc("Signed heat flux per unit area")
                )
                .bind([&model, amplitude](const std::string& target, fem::Precision value) {
                    logging::error(std::isfinite(value),
                        "HEATFLUX: value must be finite");

                    model::SurfaceRegion::Ptr region;
                    if (model._data->surface_sets.has(target)) {
                        region = model._data->surface_sets.get(target);
                    } else {
                        region = std::make_shared<model::SurfaceRegion>("INTERNAL");
                        region->add(model.compiled_surface_id(target));
                    }

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!amplitude->empty()) {
                        logging::error(model._data->amplitudes.has(*amplitude),
                            "HEATFLUX: amplitude ", *amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(*amplitude);
                    }

                    auto load = std::make_shared<bc::HeatFlux>();
                    load->region_    = std::move(region);
                    load->heat_flux_ = value;
                    load->amplitude_ = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
