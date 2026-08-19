/**
 * @file register_boundary.inl
 * @brief Registers Abaqus displacement and rotation constraints.
 *
 * Abaqus `BOUNDARY` rows are translated into six-component FEMaster supports
 * for individual nodes or compiled node sets. The command supports zero and
 * prescribed structural DOFs, applies any nodal `TRANSFORM` orientation and
 * stores the constraints in the step support collector.
 *
 * Procedure-dependent validation limits nonzero values and amplitudes to the
 * static cases that FEMaster can represent without boundary-condition history.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <limits>
#include <memory>
#include <string>

#include "../reference.h"
#include "../parser_abq.h"
#include "../../../bc/support.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_boundary(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "STEP"}));
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("DISPLACEMENT").allowed({"DISPLACEMENT"})
                .key("AMPLITUDE").optional()
        );
        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(!state.step_active || parser.active_loadcase(),
                "BOUNDARY: inside STEP must appear after a supported procedure");

            parser.model()._data->supp_cols.activate("__ABQ_STEP_SUPPORTS");
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "BOUNDARY: amplitude ", *amplitude, " does not exist");
            logging::error(state.step_active || amplitude->empty(),
                "BOUNDARY: AMPLITUDE outside STEP is not supported");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .one<int>().name("FIRST_DOF")
                    .one<int>().name("LAST_DOF").on_missing(-1).on_empty(-1)
                    .one<Precision>().name("MAGNITUDE").on_missing(Precision{0}).on_empty(Precision{0})
                )
                .bind([&parser, amplitude](const std::string& target,
                                           int first_dof,
                                           int last_dof,
                                           Precision magnitude) {
                    if (last_dof < 0) last_dof = first_dof;
                    logging::error(first_dof >= 1 && first_dof <= 6
                                && last_dof >= first_dof && last_dof <= 6,
                        "BOUNDARY: structural DOFs must be in [1,6]");

                    if (parser.abaqus_state().step_active && magnitude != Precision(0)) {
                        const std::string procedure = parser.active_loadcase()->type_name();
                        logging::error(procedure == "LINEARSTATIC"
                                    || procedure == "NONLINEARSTATIC",
                            "BOUNDARY: nonzero values are supported only for static procedures");
                        if (!amplitude->empty()) {
                            logging::error(procedure == "LINEARSTATIC",
                                "BOUNDARY: nonzero AMPLITUDE is supported only for linear static procedures");
                            magnitude *= parser.model()._data->amplitudes.get(*amplitude)->evaluate(
                                parser.abaqus_state().step_period);
                        }
                    }

                    Vec6 values;
                    values.setConstant(std::numeric_limits<Precision>::quiet_NaN());
                    for (int dof = first_dof; dof <= last_dof; ++dof) values[dof - 1] = magnitude;

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    const auto add_node = [&](ID node_id) {
                        cos::CoordinateSystem::Ptr orientation = nullptr;
                        const auto transform = state.node_transforms.find(node_id);
                        if (transform != state.node_transforms.end()) {
                            logging::error(model._data->coordinate_systems.has(transform->second),
                                "BOUNDARY: coordinate system ", transform->second, " does not exist");
                            orientation = model._data->coordinate_systems.get(transform->second);
                        }

                        auto region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(node_id);
                        model.add_support(bc::Support{std::move(region), values, std::move(orientation)});
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
