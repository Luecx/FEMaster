/**
 * @file register_cload.cpp
 * @brief Translate Abaqus nodal loads using the shared CLOAD row format.
 *
 * Standard Abaqus DOF rows remain supported. FEMaster additionally accepts
 * six-component vector rows and named model-level CLOAD definitions.
 */
#include "../commands/register_functions.h"
#include "../commands/cload_common.h"
#include "../../dsl/registry.h"

#include <array>
#include <memory>
#include <string>

#include "../parser_abq.h"
#include "../../../bc/neumann/load_c.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

void register_cload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    namespace dsl = fem::io::dsl;
    registry.command("CLOAD", [&](dsl::Command& command) {
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "STEP"}));
        auto amplitude = std::make_shared<std::string>();
        auto in_step   = std::make_shared<bool>(false);

        command.keyword(
            dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").alternative("LOADCOLLECTOR").alternative("NAME").optional()
                .key("AMPLITUDE").optional()
                .flag("FOLLOWER")
                .flag("REAL")
                .flag("IMAGINARY")
        );
        command.on_enter([&parser, amplitude, in_step](const dsl::ParentInfo& parent,
                                                      const dsl::Keys& keys) {
            *in_step = parent.command == "STEP";
            if (*in_step) {
                auto* loadcase = parser.active_loadcase();
                logging::error(parser.abaqus_state().step_active && loadcase != nullptr,
                    "CLOAD: must appear after a supported procedure inside STEP");
                logging::error(loadcase->type_name() != "EIGENFREQ",
                    "CLOAD: not supported in a FREQUENCY step");
            }
            logging::error(!keys.has("FOLLOWER"), "CLOAD: FOLLOWER is not supported");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "CLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"), "CLOAD: IMAGINARY is not supported");
            *amplitude = keys.raw("AMPLITUDE");
            if (!*in_step && !amplitude->empty()) {
                logging::error(parser.model()._data->amplitudes.has(*amplitude),
                    "CLOAD: unknown amplitude ", *amplitude);
            }
            parser.activate_cload_collector(keys.raw("LOAD_COLLECTOR"), *in_step,
                                            "__ABQ_STEP_LOADS");
        });

        // Other Abaqus load commands use the step's default active collector.
        // A named CLOAD must not redirect subsequent DLOAD/DSLOAD entries.
        command.on_exit([&parser, in_step](const dsl::Keys&) {
            if (*in_step) parser.model()._data->load_cols.activate("__ABQ_STEP_LOADS");
        });

        // Read one record at a time; format can change between rows.
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .fixed<std::string, 6>().name("COMPONENTS")
                        .desc("DOF,magnitude or Fx,Fy,Fz[,Mx,My,Mz]")
                        .on_missing(std::string{commands::cload_common::missing_token})
                        .on_empty(std::string{commands::cload_common::empty_token})
                )
                .bind([&parser, amplitude, in_step](const std::string& target,
                                                   const std::array<std::string, 6>& components) {
                    auto values = commands::cload_common::parse(components);
                    auto& model = parser.model();
                    std::string resolved_amplitude = *amplitude;
                    Precision scale = Precision(1);
                    if (*in_step) {
                        auto resolved = parser.resolve_load_amplitude(*amplitude);
                        scale = resolved.first;
                        resolved_amplitude = std::move(resolved.second);
                    }
                    values *= scale;
                    if (values.isZero()) return;

                    bc::Amplitude::Ptr load_amplitude = nullptr;
                    if (!resolved_amplitude.empty()) {
                        logging::error(model._data->amplitudes.has(resolved_amplitude),
                            "CLOAD: amplitude ", resolved_amplitude, " does not exist");
                        load_amplitude = model._data->amplitudes.get(resolved_amplitude);
                    }

                    auto& state = parser.abaqus_state();
                    const auto add_node = [&](ID node_id) {
                        cos::CoordinateSystem::Ptr orientation = nullptr;
                        const auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            logging::error(model._data->coordinate_systems.has(transform->second),
                                "CLOAD: coordinate system ", transform->second, " does not exist");
                            orientation = model._data->coordinate_systems.get(transform->second);
                        }
                        auto region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(node_id);
                        commands::cload_common::add(model, std::move(region), values,
                                                   std::move(orientation), load_amplitude);
                    };

                    if (model._data->node_sets.has(target)) {
                        for (const ID node_id : *model._data->node_sets.get(target)) add_node(node_id);
                    } else {
                        add_node(model.compiled_node_id(target));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
