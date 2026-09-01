/**
 * @file register_cload.cpp
 * @brief Registers the CLOAD input command.
 *
 * Each `CLOAD` data row creates one `bc::CLoad` for a compiled node reference or
 * node set. The six generalized components are ordered as
 * `[Fx, Fy, Fz, Mx, My, Mz]`; missing or empty components default to zero.
 *
 * The command selects a named load collector and may resolve a coordinate system
 * and amplitude shared by all of its data rows. Instance-qualified scalar node
 * references are mapped into compiled assembly identifiers. Coordinate-system
 * transformation and amplitude evaluation remain the responsibility of
 * `bc::CLoad::apply()` during load assembly.
 *
 * @see bc::CLoad
 * @see model::Model::add_load
 * @see model::Model::compiled_node_id
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <memory>
#include <string>

#include "../../../bc/neumann/load_c.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers concentrated nodal force and moment assignments.
 *
 * Command entry resolves the optional coordinate system and amplitude once and
 * activates the requested load collector. Every subsequent data row resolves
 * its node target, creates one six-component `bc::CLoad` and transfers it to
 * that collector through `Model::add_load()`.
 *
 * @param registry Parser registry receiving the command definition.
 * @param model Compiled model providing node regions, shared definitions and
 *              load collectors.
 */
void register_cload(dsl::Registry& registry, model::Model& model) {
    registry.command("CLOAD", [&](dsl::Command& command) {
        // Concentrated loads operate on compiled assembly nodes at root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Expose the component order and optional modifiers through the registry documentation
        command.doc(
            "Create concentrated nodal loads ordered as Fx, Fy, Fz, Mx, My, Mz. "
            "Optional orientation and amplitude references apply to every data row."
        );

        // Retain resolved command-wide modifiers across all data-row callbacks
        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr        >(nullptr);

        // Define the target collector and optional load modifiers
        command.keyword(
            dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .alternative("LOADCOLLECTOR")
                    .alternative("NAME")
                    .required()
                    .doc("Collector receiving the concentrated loads")
                .key("ORIENTATION").optional().doc("Coordinate system for the load components")
                .key("AMPLITUDE"  ).optional().doc("Amplitude scaling the complete generalized load")
        );

        // Resolve shared references and select the collector for this command occurrence
        command.on_enter([&model, orientation, amplitude](const dsl::Keys& keys) {
            orientation->reset();
            amplitude  ->reset();

            const std::string orientation_name = keys.raw("ORIENTATION");
            const std::string amplitude_name   = keys.raw("AMPLITUDE");

            if (!orientation_name.empty()) {
                logging::error(model._data->coordinate_systems.has(orientation_name),
                    "CLOAD: coordinate system ", orientation_name, " does not exist");
                *orientation = model._data->coordinate_systems.get(orientation_name);
            }

            if (!amplitude_name.empty()) {
                logging::error(model._data->amplitudes.has(amplitude_name),
                    "CLOAD: amplitude ", amplitude_name, " does not exist");
                *amplitude = model._data->amplitudes.get(amplitude_name);
            }

            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));
        });

        // Read one node target and six generalized load components per data row
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>   ().name("TARGET").desc("Compiled node set or scalar node reference")
                    .fixed<Precision, 6>().name("LOAD"  ).desc("Fx, Fy, Fz, Mx, My, Mz")
                        .on_missing(Precision{0}).on_empty(Precision{0})
                )
                .bind([&model, orientation, amplitude](const std::string& target,
                                                       const std::array<Precision, 6>& values) {
                    // Reuse a named compiled node set or create a private single-node region
                    model::NodeRegion::Ptr region = nullptr;
                    if (model._data->node_sets.has(target)) {
                        region = model._data->node_sets.get(target);
                    } else {
                        region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(model.compiled_node_id(target));
                    }

                    // Assemble the nominal load and retain its optional runtime modifiers
                    auto load = std::make_shared<bc::CLoad>();
                    load->region_      = std::move(region);
                    load->orientation_ = *orientation;
                    load->amplitude_   = *amplitude;
                    load->values_ << values[0], values[1], values[2], values[3], values[4], values[5];

                    // Transfer the completed load to the collector activated on command entry
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
