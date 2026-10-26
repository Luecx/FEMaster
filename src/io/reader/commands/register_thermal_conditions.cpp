/**
 * @file register_thermal_conditions.cpp
 * @brief Registers persistent scalar thermal boundary-condition histories.
 *
 * TEMPERATURE defines structural target-temperature fields. HEATFLUX and
 * CONVECTION remain thermal condition histories for heat-transfer analyses.
 * OP=NEW on TEMPERATURE resets only the current step's target field.
 *
 * Thermal boundary conditions are distinct from structural predefined fields.
 *
 * @see bc::HeatFlux
 * @see bc::Convection
 * @see bc::ConditionState
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../../bc/thermal/convection.h"
#include "../../../bc/thermal/heat_flux.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

/**
 * Registers the three independent thermal input-history families.
 *
 * TEMPERATURE writes nodal structural target fields. HEATFLUX and CONVECTION
 * remain separate thermal condition histories for future heat-transfer steps.
 *
 * @param registry Registry receiving the thermal command grammar.
 * @param parser Parser exposing the authoritative model state.
 */
void register_thermal_conditions(fem::io::dsl::Registry& registry, Parser& parser) {
    namespace dsl = fem::io::dsl;

    // Predefined temperatures are structural field values, not constraints on
    // thermal DOFs. Heat-transfer Dirichlet conditions remain separate.
    registry.command("TEMPERATURE", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"LOADCASE", "STATIC"}));
        command.doc("Set absolute nodal target temperatures for a structural static step.");
        command.keyword(dsl::KeywordSpec::make()
            .key("OP").optional("MOD").allowed({"MOD", "NEW"})
        );
        command.on_enter([&parser](const dsl::Keys& keys) {
            auto* active = parser.active_loadcase();
            logging::error(active != nullptr &&
                           (active->type_name() == "LINEARSTATIC" ||
                            active->type_name() == "NONLINEARSTATIC"),
                "TEMPERATURE: requires a structural static analysis");
            if (keys.raw("OP") == "NEW") {
                auto& data = *parser.model()._data;
                const auto& base = active->resume
                    ? data.temperature : data.temperature_reference;
                data.temperature_target = base
                    ? std::make_shared<model::Field>(*base) : nullptr;
            }
        });
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled node set or scalar node reference")
                    .one<Precision>().name("TEMPERATURE").desc("Absolute target temperature")
                )
                .bind([&parser](const std::string& target, Precision value) {
                    logging::error(std::isfinite(value),
                        "TEMPERATURE: prescribed value must be finite");
                    auto& data = *parser.model()._data;
                    if (!data.temperature_target) {
                        data.temperature_target = std::make_shared<model::Field>(
                            "TEMPERATURE_TARGET", model::FieldDomain::NODE,
                            data.field_rows(model::FieldDomain::NODE), 1);
                        data.temperature_target->fill_nan();
                    }
                    const auto region = parser.model().resolve_node_region(target);
                    for (const ID id : region->data()) {
                        logging::error(id >= 0 &&
                            static_cast<Index>(id) < data.temperature_target->rows,
                            "TEMPERATURE: node index out of bounds: ", id);
                        (*data.temperature_target)(static_cast<Index>(id), 0) = value;
                    }
                })
            )
        );

        // Empty OP=NEW blocks reset this step's target to its starting field.
        command.variant(dsl::Variant::make()
            .when(dsl::Condition::key_equals("OP", "NEW"))
        );
    });

    // Prescribed heat density supplies the scalar thermal right-hand side
    registry.command("HEATFLUX", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Prescribe surface heat input per unit reference area.");
        auto amplitude = std::make_shared<bc::Amplitude::Ptr>(nullptr);
        command.keyword(dsl::KeywordSpec::make()
            .key("OP").optional("MOD").allowed({"MOD", "NEW"})
            .key("AMPLITUDE").optional()
        );
        command.on_enter([&parser, amplitude](const dsl::Keys& keys) {
            // Resolve the optional scalar history before modifying active definitions
            auto& data = *parser.model()._data;
            const std::string name = keys.raw("AMPLITUDE");
            logging::error(name.empty() || data.amplitudes.has(name),
                "HEATFLUX: amplitude ", name, " does not exist");
            amplitude->reset();
            if (!name.empty()) {
                *amplitude = data.amplitudes.get(name);
            }
            if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::HEAT_FLUX);
            }
        });
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled surface set or scalar surface reference")
                    .one<Precision>().name("HEAT_FLUX").desc("Heat density positive into the thermal balance")
                )
                .bind([&parser, amplitude](const std::string& target, Precision value) {
                    // Keep nominal heat input and shared amplitude for later assembly
                    auto heat_flux = std::make_shared<bc::HeatFlux>();
                    heat_flux->region_    = parser.model().resolve_surface_region(target);
                    heat_flux->heat_flux_ = value;
                    heat_flux->amplitude_ = *amplitude;
                    parser.modify_conditions(bc::HEAT_FLUX, "SURFACE:" + target, {std::move(heat_flux)});
                })
            )
        );

        // Empty OP=NEW blocks are valid resets; on_enter handles the history update.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("OP", "NEW"))
        );
    });

    // Newton cooling contributes both the ambient source and the film operator
    registry.command("CONVECTION", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Define surface convection by film coefficient and ambient temperature.");
        auto amplitude = std::make_shared<bc::Amplitude::Ptr>(nullptr);
        command.keyword(dsl::KeywordSpec::make()
            .key("OP").optional("MOD").allowed({"MOD", "NEW"})
            .key("AMPLITUDE").optional()
        );
        command.on_enter([&parser, amplitude](const dsl::Keys& keys) {
            // Amplitude scales h; the prescribed ambient temperature stays absolute
            auto& data = *parser.model()._data;
            const std::string name = keys.raw("AMPLITUDE");
            logging::error(name.empty() || data.amplitudes.has(name),
                "CONVECTION: amplitude ", name, " does not exist");
            amplitude->reset();
            if (!name.empty()) {
                *amplitude = data.amplitudes.get(name);
            }
            if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::CONVECTION);
            }
        });
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled surface set or scalar surface reference")
                    .one<Precision>().name("FILM_COEFFICIENT").desc("Non-negative Newton cooling coefficient")
                    .one<Precision>().name("AMBIENT_TEMPERATURE").desc("Absolute ambient temperature")
                )
                .bind([&parser, amplitude](const std::string& target, Precision film, Precision ambient) {
                    // Replace the surface cooling prescription within its own family
                    auto convection = std::make_shared<bc::Convection>();
                    convection->region_              = parser.model().resolve_surface_region(target);
                    convection->film_coefficient_    = film;
                    convection->ambient_temperature_ = ambient;
                    convection->amplitude_            = *amplitude;
                    parser.modify_conditions(bc::CONVECTION, "SURFACE:" + target, {std::move(convection)});
                })
            )
        );

        // Empty OP=NEW blocks are valid resets; on_enter handles the history update.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("OP", "NEW"))
        );
    });
}

} // namespace fem::io::reader::commands
