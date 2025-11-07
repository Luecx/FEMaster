// register_cload.inl — DSL registration for *CLOAD
#pragma once
/**
 * @file register_cload.inl
 * @brief Register *CLOAD command.
 *
 * Defines the user-facing syntax for prescribing concentrated nodal loads.
 * Each record accepts force (Fx,Fy,Fz) and moment (Mx,My,Mz) components and
 * can optionally reference an orientation coordinate system. When an
 * orientation is supplied the components are interpreted in that local frame
 * and rotated into the global basis prior to assembly, enabling cylindrical or
 * arbitrarily rotated load directions without manual projection.
 */

#include <array>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_cload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CLOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Apply concentrated nodal forces and moments. Each data line supplies six components "
            "(Fx,Fy,Fz,Mx,My,Mz) for a target node or node set. Values accumulate when multiple "
            "records address the same entity. Optional ORIENTATION lets you specify the components "
            "in a local coordinate system, while AMPLITUDE references a time history used to scale "
            "the load during transient analyses."
        );

        auto orientation = std::make_shared<std::string>();
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                    .doc("Name of the active load collector that receives the contributions.")
                .key("ORIENTATION").optional()
                    .doc("Optional coordinate system used to interpret the six load components.")
                .key("AMPLITUDE").optional()
                    .doc("Optional time amplitude that scales the load components.")
        );

        command.on_enter([orientation, amplitude, &model](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(collector);
        });

        command.variant(
            fem::dsl::Variant::make()
                .doc("Each data line defines a target node set/id followed by Fx, Fy, Fz, Mx, My, Mz (local or global depending on ORIENTATION).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(1)) // one or more lines
                        .pattern(
                            fem::dsl::Pattern::make()
                                .one<std::string>()
                                    .name("TARGET")
                                    .desc("Node set name or single node id (integer).")
                                .fixed<fem::Precision, 6>()
                                    .name("LOAD")
                                    .desc("Fx, Fy, Fz, Mx, My, Mz. Uses the ORIENTATION basis when provided, otherwise global axes.")
                                    .on_missing(fem::Precision{0})
                                    .on_empty  (fem::Precision{0})
                        )
                        .bind([&model, orientation, amplitude](const std::string& target,
                                       const std::array<fem::Precision, 6>& values) {
                            fem::Vec6 load;
                            load << values[0], values[1], values[2],
                                    values[3], values[4], values[5];

                            // node set target?
                            if (model._data->node_sets.has(target)) {
                                model.add_cload(target, load, *orientation, *amplitude);
                                return;
                            }

                            // single node id?
                            try {
                                const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                                model.add_cload(id, load, *orientation, *amplitude);
                            } catch (const std::exception&) {
                                throw std::runtime_error(
                                    "CLOAD target '" + target +
                                    "' is neither an existing node set nor a valid node id."
                                );
                            }
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
