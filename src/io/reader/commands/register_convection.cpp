/**
 * @file register_convection.cpp
 * @brief Registers linear thermal convection boundary conditions.
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "../../../bc/robin/convection.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

void register_convection(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONVECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Apply q = h (T_inf - T) Robin convection to compiled surfaces.");

        auto amplitude = std::make_shared<std::string>();
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .alternative("LOAD COLLECTOR")
                    .alternative("NAME")
                    .required()
                    .doc("Load collector receiving the convection conditions")
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
                    .one<fem::Precision>().name("FILM_COEFFICIENT").desc("Convection film coefficient h")
                    .one<fem::Precision>().name("AMBIENT_TEMPERATURE").desc("Prescribed ambient temperature T_inf")
                )
                .bind([&model, amplitude](const std::string& target,
                                          fem::Precision film_coefficient,
                                          fem::Precision ambient_temperature) {
                    logging::error(std::isfinite(film_coefficient) && film_coefficient >= Precision(0),
                        "CONVECTION: film coefficient must be finite and non-negative");
                    logging::error(std::isfinite(ambient_temperature),
                        "CONVECTION: ambient temperature must be finite");

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
                            "CONVECTION: amplitude ", *amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(*amplitude);
                    }

                    auto load = std::make_shared<bc::Convection>();
                    load->region_              = std::move(region);
                    load->film_coefficient_    = film_coefficient;
                    load->ambient_temperature_ = ambient_temperature;
                    load->amplitude_           = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
