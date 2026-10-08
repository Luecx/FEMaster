/**
 * @file register_vload.cpp
 * @brief Registers FEMaster volume force densities.
 *
 * Each row resolves a compiled element region through model::Model and creates
 * a bc::VLoad. Named definitions populate only their explicit collector;
 * unnamed analysis commands apply VLOAD only during their active step. The DSL owns
 * MOD/NEW replacement and LOADS selection. The load implementation owns volume
 * integration, coordinate transformation and amplitude evaluation. Nodal
 * TRANSFORM does not define the basis of distributed loads.
 *
 * @see bc::VLoad
 * @see Parser::activate_load_collector
 * @see model::Model::resolve_element_region
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../parser.h"
#include "../../../bc/structural/load_v.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <array>
#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers volume force densities and their collector context.
 *
 * Entry resolves the optional orientation and amplitude once for all data rows and
 * selects named storage when NAME is supplied. Unnamed analysis rows update
 * temporary VLOADs; OP=NEW also retires inherited direct definitions in the
 * family before reading the block. NAME is required outside an analysis.
 * Each target is a named compiled element set or a scalar, possibly instance-
 * qualified reference. Missing or empty load components retain their zero default.
 * Integration and conversion into consistent nodal forces remain in bc::VLoad.
 *
 * @param registry Registry receiving the command grammar and callbacks.
 * @param parser Parser supplying the compiled model and active analysis.
 */
void register_vload(dsl::Registry& registry, Parser& parser) {
    registry.command("VLOAD", [&](dsl::Command& command) {
        // Distributed loads operate on compiled regions in model or analysis scope.
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Create FEMaster volume force densities; Abaqus-compatible BX/BY/BZ are available through DLOAD.");

        // Resolved modifiers are reset on entry and shared by this occurrence's rows.
        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr>(nullptr);

        // An explicit NAME selects reusable definition storage. Unnamed analysis
        // rows update direct history and never inherit a prior insertion target.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME").alternative("LOAD_COLLECTOR").alternative("LOADCOLLECTOR")
                    .optional().doc("Collector name; optional inside LOADCASE")
                .key("ORIENTATION").optional().doc("Coordinate system for the load components")
                .key("AMPLITUDE").optional().doc("Amplitude scaling the load")
        );
        command.on_enter([&parser, orientation, amplitude, collector](const dsl::Keys& keys) {
            auto& model = parser.model();

            orientation->reset();
            amplitude->reset();

            // Resolve reusable model definitions before changing the active collector.
            const std::string orientation_name = keys.raw("ORIENTATION");
            const std::string amplitude_name   = keys.raw("AMPLITUDE");

            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "VLOAD: coordinate system ", orientation_name, " does not exist");
            logging::error(amplitude_name.empty() || model._data->amplitudes.has(amplitude_name),
                "VLOAD: amplitude ", amplitude_name, " does not exist");

            // Resolve validated references; omitted modifiers retain their null handles.
            if (!orientation_name.empty()) {
                *orientation = model._data->coordinate_systems.get(orientation_name);
            }
            if (!amplitude_name.empty()) {
                *amplitude = model._data->amplitudes.get(amplitude_name);
            }
            // NAME denotes a named definition, while an unnamed row in an
            // analysis changes only the current model history.
            const std::string collector_name = keys.raw("NAME");
            logging::error(!collector_name.empty() || parser.active_loadcase() != nullptr,
                "VLOAD: definitions outside an analysis require NAME");
            collector->reset();
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            } else if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::VLOAD);
            }
        });

        // Keep the body-force-density field in its element region. Nodal
        // TRANSFORM assignments apply to nodal quantities and do not alter it.
        command.data(
            dsl::Pattern::make()
                .one<std::string   >("TARGET", "Compiled element set or scalar reference")
                .fixed<Precision, 3>("LOAD", "Three load components")
                    .defaults(Precision{0}),
            [&parser, orientation, amplitude, collector](
                const std::string&              target,
                const std::array<Precision, 3>& values
            ) {
                auto& model = parser.model();

                // Resolve the target through the model and retain nominal values
                // and shared modifiers for integration during load assembly.
                auto load = std::make_shared<bc::VLoad>();
                load->region_      = model.resolve_element_region(target);
                load->values_      = Vec3{values[0], values[1], values[2]};
                load->orientation_ = *orientation;
                load->amplitude_   = *amplitude;

                // Named definitions and direct history are mutually exclusive targets
                if (*collector) {
                    (*collector)->add(std::move(load));
                } else {
                    parser.select_collector_condition(bc::VLOAD, std::move(load));
                }
            }
        );
    });
}

} // namespace fem::io::reader::commands
