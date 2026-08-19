/**
 * @file register_cload.inl
 * @brief Registers Abaqus concentrated nodal loads.
 *
 * Each Abaqus `CLOAD` row maps one structural DOF and magnitude onto a
 * six-component FEMaster concentrated load for a node or node set. Nodal
 * `TRANSFORM` assignments define the optional load basis, and supported
 * amplitude semantics are resolved against the active analysis procedure.
 *
 * Real-valued non-follower loading is stored in the dedicated step load
 * collector. Unsupported complex or follower variants are rejected explicitly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>
#include <string>

#include "../reference.h"
#include "../parser_abq.h"
#include "../../../bc/load_c.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_cload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("CLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional()
                .flag("FOLLOWER")
                .flag("REAL")
                .flag("IMAGINARY")
        );
        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto* loadcase = parser.active_loadcase();
            logging::error(parser.abaqus_state().step_active && loadcase != nullptr,
                "CLOAD: must appear after a supported procedure inside STEP");
            logging::error(loadcase->type_name() != "EIGENFREQ",
                "CLOAD: not supported in a FREQUENCY step");
            logging::error(!keys.has("FOLLOWER"),
                "CLOAD: FOLLOWER is not supported");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "CLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "CLOAD: IMAGINARY is not supported");
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<int>().name("DOF")
                    .one<Precision>().name("MAGNITUDE")
                )
                .bind([&parser, amplitude](const std::string& target, int dof, Precision magnitude) {
                    logging::error(dof >= 1 && dof <= 6,
                        "CLOAD: DOF must be in [1,6]");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    magnitude *= scale;
                    if (magnitude == Precision(0)) return;

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    Vec6 values = Vec6::Zero();
                    values[dof - 1] = magnitude;

                    const auto add_node = [&](ID node_id) {
                        cos::CoordinateSystem::Ptr orientation = nullptr;
                        const auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            logging::error(model._data->coordinate_systems.has(transform->second),
                                "CLOAD: coordinate system ", transform->second, " does not exist");
                            orientation = model._data->coordinate_systems.get(transform->second);
                        }

                        bc::Amplitude::Ptr load_amplitude = nullptr;
                        if (!resolved_amplitude.empty()) {
                            logging::error(model._data->amplitudes.has(resolved_amplitude),
                                "CLOAD: amplitude ", resolved_amplitude, " does not exist");
                            load_amplitude = model._data->amplitudes.get(resolved_amplitude);
                        }

                        auto region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(node_id);
                        auto load = std::make_shared<bc::CLoad>();
                        load->region_      = std::move(region);
                        load->values_      = values;
                        load->orientation_ = std::move(orientation);
                        load->amplitude_   = std::move(load_amplitude);
                        model.add_load(std::move(load));
                    };

                    if (model._data->node_sets.has(target)) {
                        for (const ID node_id : *model._data->node_sets.get(target)) add_node(node_id);
                    } else {
                        add_node(io::reader::compiled_node_id(model, target));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
