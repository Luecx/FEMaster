/**
 * @file register_boundary.cpp
 * @brief Registers Abaqus mechanical and scalar-temperature constraints.
 *
 * Abaqus `BOUNDARY` rows are translated into shared FEMaster support
 * definitions for individual nodes or compiled node sets. Every prescribed
 * generalized component becomes one condition object, so later history updates
 * can replace one DOF without affecting independent support components. Nodal
 * `TRANSFORM` assignments provide the optional local basis.
 *
 * Mechanical DOFs 1-6 create SUPPORT definitions, while Abaqus temperature
 * DOF 11 creates independent TEMPERATURE constraints for heat-transfer steps.
 * OP=MOD replaces matching targets; OP=NEW resets the current procedure's
 * boundary-condition family. Initial BOUNDARY definitions seed the matching
 * persistent history independently of the selected analysis procedure.
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
#include "../../../bc/thermal/temperature.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers initial and STEP-local prescribed mechanical or thermal components.
 *
 * Structural DOFs retain the existing support splitting and nodal TRANSFORM
 * semantics. DOF 11 instead prescribes an absolute nodal scalar temperature;
 * it is mapped to thermal component zero when the condition is assembled.
 * The two independent condition families preserve their own history updates.
 *
 * @param registry Abaqus registry receiving the BOUNDARY command.
 * @param parser Reader supplying compiled targets and current procedure controls.
 */
void register_boundary(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("BOUNDARY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE",
            "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS", "HEATTRANSFER"}));
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
                // A thermal step replaces temperature constraints, not structural supports.
                // Outside a step, the keyword may have initialized either family.
                const auto* active = parser.active_loadcase();
                const bool thermal = active != nullptr && active->type_name() == "STEADYHEATTRANSFER";
                if (!thermal) parser.clear_conditions(bc::SUPPORT);
                if (thermal || active == nullptr) parser.clear_conditions(bc::TEMPERATURE);
            }
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
                    logging::error((first_dof >= 1 && first_dof <= 6
                                    && last_dof >= first_dof && last_dof <= 6)
                                || (first_dof == 11 && last_dof == 11),
                        "BOUNDARY: valid components are structural DOFs 1-6 or temperature DOF 11");
                    logging::error(std::isfinite(magnitude),
                        "BOUNDARY: prescribed value must be finite");

                    const auto* active = parser.active_loadcase();
                    const bool thermal = active != nullptr && active->type_name() == "STEADYHEATTRANSFER";
                    logging::error(active == nullptr || (thermal == (first_dof == 11)),
                        "BOUNDARY: temperature DOF 11 is available only in HEAT TRANSFER, "
                        "and structural DOFs 1-6 only in structural steps");

                    if (parser.step_state().step_active && magnitude != Precision(0)) {
                        const std::string procedure = parser.active_loadcase()->type_name();
                        logging::error(procedure == "LINEARSTATIC"
                                    || procedure == "NONLINEARSTATIC"
                                    || procedure == "STEADYHEATTRANSFER",
                            "BOUNDARY: nonzero values are supported only for static procedures");
                        if (!amplitude->empty()) {
                            logging::error(procedure == "LINEARSTATIC"
                                        || procedure == "STEADYHEATTRANSFER",
                                "BOUNDARY: nonzero AMPLITUDE requires a linear procedure");
                            magnitude *= parser.model()._data->amplitudes.get(*amplitude)->evaluate(
                                parser.step_state().step_period);
                        }
                    }

                    auto& model = parser.model();
                    const std::string identifier =
                        (model._data->node_sets.has(target) ? "NSET:" : "NODE:") + target;

                    // Abaqus temperature DOF 11 maps to local scalar thermal DOF 0.
                    // Keep it separate from generalized mechanical support components.
                    if (first_dof == 11) {
                        auto temperature = std::make_shared<bc::Temperature>(
                            model.resolve_node_region(target), magnitude);
                        parser.modify_conditions(bc::TEMPERATURE, identifier, {std::move(temperature)});
                        return;
                    }

                    Vec6 values;
                    values.setConstant(std::numeric_limits<Precision>::quiet_NaN());
                    for (int dof = first_dof; dof <= last_dof; ++dof) values[dof - 1] = magnitude;

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
                })
            )
        );

        // Abaqus/CAE can emit OP=NEW without data rows to clear prior constraints.
        // Keep all other BOUNDARY occurrences subject to the normal row requirement.
        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("OP", "NEW"))
        );
    });
}

} // namespace fem::io::reader::commands