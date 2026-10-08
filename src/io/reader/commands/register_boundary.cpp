/**
 * @file register_boundary.cpp
 * @brief Registers Abaqus displacement and rotation constraints.
 *
 * Abaqus `BOUNDARY` rows are translated into shared FEMaster support
 * definitions for individual nodes or compiled node sets. Every prescribed
 * generalized component becomes one condition object, so later history updates
 * can replace one DOF without affecting independent support components. Nodal
 * `TRANSFORM` assignments provide the optional local basis.
 *
 * The resulting shared supports are retained in ModelData's persistent
 * SUPPORT history. OP=MOD replaces matching node/DOF definitions, while OP=NEW
 * clears the complete support family before the rows in the current block are
 * applied. Initial BOUNDARY definitions outside STEP seed the same history.
 *
 * Procedure-dependent validation of nonzero values and amplitudes remains tied
 * to the currently supported solver behavior.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../commands/register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "../parser.h"
#include "../../../bc/structural/support.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers initial and STEP-local prescribed structural components.
 *
 * Initial and analysis definitions modify the same ModelData SUPPORT history.
 * NEW clears this family when the keyword is encountered; MOD replaces only
 * matching normalized node/DOF targets. Each row is split into individual
 * prescriptions and nodal TRANSFORM bases are retained for equation assembly.
 * Supported procedures validate nonzero displacements and amplitude scaling
 * before definitions are added. END STEP reads the resulting history directly.
 *
 * @param registry Abaqus registry receiving the BOUNDARY command.
 * @param parser Reader supplying compiled targets and current procedure controls.
 */
void register_boundary(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE",
            "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("DISPLACEMENT").allowed({"DISPLACEMENT"})
                .key("AMPLITUDE").optional()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
        );
        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.step_state();
            logging::error(!state.step_active || parser.active_loadcase(),
                "BOUNDARY: inside STEP must appear after a supported procedure");

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "BOUNDARY: amplitude ", *amplitude, " does not exist");
            logging::error(state.step_active || amplitude->empty(),
                "BOUNDARY: AMPLITUDE outside STEP is not supported");

            // Initial and step-local boundaries share one persistent history.
            // NEW resets that target history immediately when the keyword is
            // encountered; MOD leaves unrelated support slots untouched.
            if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::SUPPORT);
            }
        });

        command.data(
            fem::io::dsl::Pattern::make()
                .one<std::string>("TARGET")
                .one<int>("FIRST_DOF")
                .one<int>("LAST_DOF").defaults(-1)
                .one<Precision>("MAGNITUDE").defaults(Precision{0}),
            [&parser, amplitude](const std::string& target,
                                       int first_dof,
                                       int last_dof,
                                       Precision magnitude) {
                if (last_dof < 0) last_dof = first_dof;
                logging::error(first_dof >= 1 && first_dof <= 6
                            && last_dof >= first_dof && last_dof <= 6,
                    "BOUNDARY: structural DOFs must be in [1,6]");

                if (parser.step_state().step_active && magnitude != Precision(0)) {
                    const std::string procedure = parser.active_loadcase()->type_name();
                    logging::error(procedure == "LINEARSTATIC"
                                || procedure == "NONLINEARSTATIC",
                        "BOUNDARY: nonzero values are supported only for static procedures");
                    if (!amplitude->empty()) {
                        logging::error(procedure == "LINEARSTATIC",
                            "BOUNDARY: nonzero AMPLITUDE is supported only for linear static procedures");
                        magnitude *= parser.model()._data->amplitudes.get(*amplitude)->evaluate(
                            parser.step_state().step_period);
                    }
                }

                Vec6 values;
                values.setConstant(std::numeric_limits<Precision>::quiet_NaN());
                for (int dof = first_dof; dof <= last_dof; ++dof) values[dof - 1] = magnitude;

                auto& model = parser.model();
                const std::string identifier =
                    (model._data->node_sets.has(target) ? "NSET:" : "NODE:") + target;
                std::array<std::vector<bc::Condition::Ptr>, 6> replacements;

                const auto add_node = [&](ID node_id) {
                    cos::CoordinateSystem::Ptr orientation = nullptr;
                    const auto transform = parser.node_transforms.find(node_id);
                    if (transform != parser.node_transforms.end()) {
                        logging::error(model._data->coordinate_systems.has(transform->second),
                            "BOUNDARY: coordinate system ", transform->second, " does not exist");
                        orientation = model._data->coordinate_systems.get(transform->second);
                    }

                    auto region = std::make_shared<model::NodeRegion>("INTERNAL");
                    region->add(node_id);

                    for (Dim dof = 0; dof < 6; ++dof) {
                        if (std::isnan(values[dof])) continue;

                        Vec6 component_values = Vec6::Constant(NAN);
                        component_values[dof] = values[dof];
                        replacements[dof].push_back(std::make_shared<bc::Support>(
                            region, component_values, orientation
                        ));
                    }
                };

                if (model._data->node_sets.has(target)) {
                    for (const ID node_id : *model._data->node_sets.get(target)) add_node(node_id);
                } else {
                    add_node(model.compiled_node_id(target));
                }

                // A node set is one logical source for each DOF, even when
                // individual nodes use different TRANSFORM bases.
                for (Dim dof = 0; dof < 6; ++dof) {
                    if (replacements[dof].empty()) continue;
                    parser.modify_conditions(bc::SUPPORT, identifier, std::move(replacements[dof]));
                }
            }
        );
    });
}

} // namespace fem::io::reader::commands