/**
 * @file register_thermal_conditions.cpp
 * @brief Registers persistent scalar thermal boundary-condition histories.
 *
 * TEMPERATURE prescribes nodal primary variables, HEATFLUX prescribes heat input
 * per reference surface area, and CONVECTION defines Newton cooling. The native
 * and shared STEP grammar use the same ModelData::conditions as structural
 * commands. OP=MOD replaces matching targets within one input family; OP=NEW
 * retires only inherited definitions from that family before the block is read.
 *
 * Initial thermal definitions may appear outside an analysis and remain active
 * until later input modifies them. Conditions are not copied into collectors at
 * analysis completion. Named thermal grouping remains available through the
 * Model registration API independently of these direct commands.
 *
 * @see bc::Temperature
 * @see bc::HeatFlux
 * @see bc::Convection
 * @see bc::ConditionManager
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../../bc/thermal/convection.h"
#include "../../../bc/thermal/heat_flux.h"
#include "../../../bc/thermal/temperature.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

/**
 * Registers the three independent thermal input-history families.
 *
 * Each keyword directly changes the current model state. No name or private
 * collector is required. Temperatures resolve node regions and are stored as
 * absolute scalar prescriptions. Heat flux and convection resolve compiled
 * surface regions; amplitudes scale heat density or film coefficient during
 * assembly, while ambient temperature is a constant replacement value.
 *
 * @param registry Registry receiving the thermal command grammar.
 * @param parser Parser exposing the authoritative model state.
 */
void register_thermal_conditions(fem::io::dsl::Registry& registry, Parser& parser) {
    namespace dsl = fem::io::dsl;

    // Prescribed temperature changes only the scalar primary-variable history
    registry.command("TEMPERATURE", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Prescribe absolute nodal temperatures in the current model history.");
        command.keyword(dsl::KeywordSpec::make()
            .key("OP").optional("MOD").allowed({"MOD", "NEW"})
        );
        command.on_enter([&parser](const dsl::Keys& keys) {
            // NEW affects only temperature prescriptions, including initial ones
            if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::TEMPERATURE);
            }
        });
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled node set or scalar node reference")
                    .one<Precision>().name("TEMPERATURE").desc("Absolute prescribed temperature")
                )
                .bind([&parser](const std::string& target, Precision value) {
                    // Replace only the matching nodal target in temperature history
                    auto temperature = std::make_shared<bc::Temperature>(
                        parser.model().resolve_node_region(target), value
                    );
                    parser.modify_conditions(bc::TEMPERATURE,
                        (parser.model()._data->node_sets.has(target) ? "NSET:" : "NODE:") + target,
                        {std::move(temperature)});
                })
            )
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
    });
}

} // namespace fem::io::reader::commands
