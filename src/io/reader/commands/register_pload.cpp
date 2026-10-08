/**
 * @file register_pload.cpp
 * @brief Registers FEMaster scalar surface pressures.
 *
 * Each row resolves a compiled surface region through model::Model and creates
 * a bc::PLoad. Named definitions populate only their explicit collector;
 * unnamed analysis commands apply PLOAD only during their active step. The DSL owns
 * MOD/NEW replacement and LOADS selection. The load implementation owns surface
 * integration and amplitude evaluation during assembly. Nodal TRANSFORM does
 * not define the basis of distributed loads.
 *
 * @see bc::PLoad
 * @see Parser::activate_load_collector
 * @see model::Model::resolve_surface_region
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../parser.h"
#include "../../../bc/structural/load_p.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers scalar surface pressures and their collector context.
 *
 * Entry resolves the optional amplitude once for all data rows and
 * selects named storage when NAME is supplied. Unnamed analysis rows update
 * temporary PLOADs; OP=NEW also retires inherited direct definitions in the
 * family before reading the block. NAME is required outside an analysis.
 * Each target is a named compiled surface set or a scalar, possibly instance-
 * qualified reference. Missing or empty load components retain their zero default.
 * Integration and conversion into consistent nodal forces remain in bc::PLoad.
 *
 * @param registry Registry receiving the command grammar and callbacks.
 * @param parser Parser supplying the compiled model and active analysis.
 */
void register_pload(dsl::Registry& registry, Parser& parser) {
    registry.command("PLOAD", [&](dsl::Command& command) {
        // Distributed loads operate on compiled regions in model or analysis scope.
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Create FEMaster scalar surface pressures; Abaqus-compatible pressure is available through DSLOAD or DLOAD.");

        // Resolved modifiers are reset on entry and shared by this occurrence's rows.
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr>(nullptr);

        // An explicit NAME selects reusable definition storage. Unnamed analysis
        // rows update direct history and never inherit a prior insertion target.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME").alternative("LOAD_COLLECTOR").alternative("LOADCOLLECTOR")
                    .optional().doc("Collector name; optional inside LOADCASE")
                .key("AMPLITUDE").optional().doc("Amplitude scaling the load")
        );
        command.on_enter([&parser, amplitude, collector](const dsl::Keys& keys) {
            auto& model = parser.model();

            // Scalar pressure follows the deformed surface normal. Until its
            // load stiffness is implemented, reject direct nonlinear usage.
            const auto* loadcase = parser.active_loadcase();
            logging::error(loadcase == nullptr || loadcase->type_name() != "NONLINEARSTATIC",
                "PLOAD: follower pressure is not supported in nonlinear steps");

            amplitude->reset();

            // Resolve reusable model definitions before changing the active collector.
            const std::string amplitude_name = keys.raw("AMPLITUDE");
            logging::error(amplitude_name.empty() || model._data->amplitudes.has(amplitude_name),
                "PLOAD: amplitude ", amplitude_name, " does not exist");

            if (!amplitude_name.empty()) {
                *amplitude = model._data->amplitudes.get(amplitude_name);
            }
            // NAME denotes a named definition, while an unnamed row in an
            // analysis changes only the current model history.
            const std::string collector_name = keys.raw("NAME");
            logging::error(!collector_name.empty() || parser.active_loadcase() != nullptr,
                "PLOAD: definitions outside an analysis require NAME");
            collector->reset();
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            } else if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::PLOAD);
            }
        });

        // Keep distributed loads in their physical region; nodal TRANSFORM does not
        // change the prescribed field or, for pressure, its surface-normal direction.
        command.data(
            dsl::Pattern::make()
                .one<std::string>("TARGET", "Compiled surface set or scalar reference")
                .one<Precision  >("P", "Pressure positive opposite to the surface normal")
                    .defaults(Precision{0}),
            [&parser, amplitude, collector](const std::string& target, Precision value) {
                auto& model = parser.model();

                // Resolve the target through the model and retain nominal values
                // and shared modifiers for integration during load assembly.
                auto load = std::make_shared<bc::PLoad>();
                load->region_    = model.resolve_surface_region(target);
                load->pressure_  = value;
                load->amplitude_ = *amplitude;

                // Named definitions and direct history are mutually exclusive targets
                if (*collector) {
                    (*collector)->add(std::move(load));
                } else {
                    parser.select_collector_condition(bc::PLOAD, std::move(load));
                }
            }
        );
    });
}

} // namespace fem::io::reader::commands
