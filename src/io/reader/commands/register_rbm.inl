/**
 * @file register_rbm.inl
 * @brief Registers rigid-body-motion suppression.
 *
 * The parser resolves the compiled element set and appends the completed Rbm
 * object directly to ModelData. This keeps Model free of generic constraint
 * dispatch and makes the same operation equally explicit for direct C++ users.
 *
 * The target element region is resolved after compilation so rigid-body modes
 * are derived from the final assembly coordinates and active structural DOFs.
 * The concrete constraint owns equation generation and rank handling.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../constraints/types/rbm.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_rbm(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("RBM", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("ELSET").alternative("SET").optional("EALL")
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string set = keys.raw("ELSET");

            logging::error(model._data->compiled,
                "RBM: constraints require a compiled model");
            logging::error(model._data->elem_sets.has(set),
                "RBM: element set ", set, " not found");
            logging::error(model._data->elem_sets.get(set) && model._data->elem_sets.get(set)->size() > 0,
                "RBM: element set ", set, " is empty");

            model._data->rbms.emplace_back(model._data->elem_sets.get(set));
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
