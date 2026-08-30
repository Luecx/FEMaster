/**
 * @file register_convection.cpp
 * @brief Registers linear thermal convection boundary conditions.
 *
 * `CONVECTION` describes heat exchange between a compiled surface and an
 * ambient medium according to
 *
 * \f$q = h\,(T_\infty - T)\f$.
 *
 * The condition is mathematically Robin-type but is stored in FEMaster's common
 * load-side `Neumann` hierarchy. During thermal system assembly its single
 * `apply()` call contributes both the ambient source term to the right-hand side
 * and the film matrix to the left-hand side.
 *
 * Native syntax:
 *
 * @code
 * *CONVECTION, LOAD_COLLECTOR=THERMAL_LOADS
 * SURFACE_OUTER, 12.5, 293.15
 * @endcode
 *
 * Data columns are target surface, film coefficient `h`, and prescribed ambient
 * temperature `T_inf`. An optional amplitude scales the complete convection
 * condition at analysis time.
 *
 * @see bc::Convection
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

#include "../../../bc/neumann/convection.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers linear surface convection in the native DSL.
 *
 * Entering the command activates one common load collector and stores the
 * optional amplitude name. Every row resolves the target as either a compiled
 * surface set or a single compiled surface. The parsed film coefficient and
 * ambient temperature are copied directly to the resulting `bc::Convection`.
 *
 * A zero film coefficient is accepted and produces a no-op contribution; a
 * negative or non-finite coefficient is rejected because it does not represent
 * the supported passive convection law.
 *
 * @param registry DSL registry receiving the `CONVECTION` command.
 * @param model Compiled model providing surfaces, amplitudes and load storage.
 */
void register_convection(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("CONVECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Apply q = h (T_inf - T) convection to compiled surfaces.");

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
