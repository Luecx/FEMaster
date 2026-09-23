/**
 * @file register_cload.cpp
 * @brief Register the shared nodal-load data syntax in native FEMaster decks.
 *
 * A data row is either TARGET,DOF,MAGNITUDE or TARGET,Fx,Fy,Fz[,Mx,My,Mz].
 * The target collector is explicit at model/assembly level and is inferred
 * from the active load case when omitted within LOADCASE.
 */
#include "register_functions.h"
#include "cload_common.h"
#include "../parser.h"
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

void register_cload(dsl::Registry& registry, Parser& parser) {
    registry.command("CLOAD", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE"}));
        command.doc("Concentrated loads: TARGET,DOF,MAGNITUDE or TARGET,Fx,Fy,Fz[,Mx,My,Mz]. "
                    "LOAD_COLLECTOR is required outside a LOADCASE.");

        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto amplitude   = std::make_shared<bc::Amplitude::Ptr>(nullptr);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").alternative("LOADCOLLECTOR").alternative("NAME")
                    .optional().doc("Required outside LOADCASE; otherwise defaults to its inline loads")
                .key("ORIENTATION").optional().doc("Coordinate system for the load components")
                .key("AMPLITUDE").optional().doc("Amplitude scaling the complete generalized load")
        );

        command.on_enter([&parser, orientation, amplitude](const dsl::ParentInfo& parent,
                                                            const dsl::Keys& keys) {
            auto& model = parser.model();
            orientation->reset();
            amplitude->reset();

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

            parser.activate_cload_collector(keys.raw("LOAD_COLLECTOR"),
                                            parent.command == "LOADCASE");
        });

        // A single repeated segment decodes each physical row independently,
        // allowing both syntax forms even when they occur within the same block.
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Compiled node or node set")
                    .fixed<std::string, 6>().name("COMPONENTS")
                        .desc("DOF,magnitude or Fx,Fy,Fz[,Mx,My,Mz]")
                        .on_missing(std::string{cload_common::missing_token})
                        .on_empty(std::string{cload_common::empty_token})
                )
                .bind([&parser, orientation, amplitude](const std::string& target,
                                                        const std::array<std::string, 6>& components) {
                    auto& model = parser.model();
                    auto values = cload_common::parse(components);
                    model::NodeRegion::Ptr region;
                    if (model._data->node_sets.has(target)) {
                        region = model._data->node_sets.get(target);
                    } else {
                        region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(model.compiled_node_id(target));
                    }
                    cload_common::add(model, std::move(region), values, *orientation, *amplitude);
                })
            )
        );
    });
}
} // namespace fem::io::reader::commands
