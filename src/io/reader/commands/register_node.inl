/**
 * @file register_node.inl
 * @brief Registers part-local and unqualified assembly nodes.
 *
 * Node rows map sparse user identifiers to global Cartesian coordinates in the
 * currently active semantic Part. Root-level definitions are assigned to the
 * default part, while assembly-scope callbacks are coordinated separately
 * during the post-compile parser pass.
 *
 * This command retains the original identifiers and does not allocate dense
 * solver rows; that transformation belongs to `Model::compile()`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <memory>

#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_node(fem::io::dsl::Registry& registry,
                          model::Model& model,
                          std::shared_ptr<bool> assembly_scope) {
    registry.command("NODE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "PART", "ASSEMBLY"}));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").optional("NALL")
        );
        command.on_enter([&model, assembly_scope](const fem::io::dsl::Keys& keys) {
            logging::error(model._data != nullptr && !model._data->compiled,
                "NODE: nodes cannot be added after compile()");
            if (*assembly_scope) {
                model._data->parts.activate(model::Model::DEFAULT_PART_NAME);
            }
            const auto part = model._data->parts.get();
            logging::error(part != nullptr,
                "NODE: no active part is available");
            part->node_sets.activate(keys.raw("NSET"));
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<ID>().name("ID")
                    .one<Precision>().name("X").on_empty(Precision{0}).on_missing(Precision{0})
                    .one<Precision>().name("Y").on_empty(Precision{0}).on_missing(Precision{0})
                    .one<Precision>().name("Z").on_empty(Precision{0}).on_missing(Precision{0})
                )
                .bind([&model](ID id, Precision x, Precision y, Precision z) {
                    model.set_node(id, x, y, z);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
