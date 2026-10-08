/**
 * @file register_support.cpp
 * @brief Registers nodal supports.
 *
 * A `SUPPORT` row targets a node or node set and supplies six generalized
 * constraint values. Finite entries prescribe translations or rotations while
 * omitted components remain unconstrained; an optional coordinate system
 * defines the basis in which the values are expressed.
 *
 * Each finite generalized component is normalized to one shared `bc::Support`
 * definition and stored in named groups or temporary step activations.
 * Nodal TRANSFORM assignments override the command-wide ORIENTATION and split
 * node regions by their effective local bases before support creation.
 * The split preserves the assembled equations while making later condition-history replacement granular
 * at the generalized-DOF level. Actual constraint equations are assembled later
 * by the active load case.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include "../parser.h"
#include "../util.h"
#include "../../../bc/structural/support.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers reusable support groups and temporary native prescriptions.
 *
 * A non-empty NAME/SUPPORT_COLLECTOR identifies named definition storage. Without
 * a name, SUPPORT requires an active analysis and changes ModelData::conditions;
 * OP=NEW releases inherited SUPPORT history; unnamed native rows are step-local.
 * Each finite vector component becomes a separate condition so independent
 * prescribed translations and rotations remain replaceable. Equation expansion
 * and coordinate projection take place later during Model constraint assembly.
 *
 * @param registry Native command grammar receiving the registration.
 * @param parser Parser exposing the model and the active analysis scope.
 */
void register_support(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("SUPPORT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));

        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);
        auto collector   = std::make_shared<bc::SupportCollector::Ptr>(nullptr);
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("SUPPORT_COLLECTOR").optional().alternative("SUPPORT COLLECTOR").alternative("NAME")
                .key("ORIENTATION").optional()
        );
        command.on_enter([&parser, orientation, collector](const fem::io::dsl::Keys& keys) {
            auto& model = parser.model();
            orientation->reset();
            const std::string orientation_name = keys.raw("ORIENTATION");
            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "SUPPORT: coordinate system ", orientation_name, " does not exist");
            if (!orientation_name.empty()) {
                *orientation = model._data->coordinate_systems.get(orientation_name);
            }

            // Reset only inherited direct history for OP=NEW; named groups remain definitions
            const std::string name = keys.raw("SUPPORT_COLLECTOR");
            logging::error(!name.empty() || parser.active_loadcase() != nullptr,
                "SUPPORT: definitions outside an analysis require NAME or SUPPORT_COLLECTOR");
            collector->reset();
            if (!name.empty()) {
                *collector = model._data->supp_cols.activate(name);
            } else if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::SUPPORT);
            }
        });

        command.data(
            fem::io::dsl::Pattern::make()
                .one<std::string>("TARGET")
                .fixed<fem::Precision, 6>("DOF")
                    .defaults(std::numeric_limits<fem::Precision>::quiet_NaN()),
            [&parser, orientation, collector](
                const std::string&                  target,
                const std::array<fem::Precision, 6>& values
            ) {
                auto& model = parser.model();
                auto region = model.resolve_node_region(target);
                // Nodal TRANSFORM takes precedence over ORIENTATION.
                // Unnamed native supports apply only in this analysis.
                for (auto& [support_region, support_orientation] :
                     group_by_orientation(parser, std::move(region), *orientation)) {
                    for (Dim dof = 0; dof < 6; ++dof) {
                        const Precision value = values[static_cast<std::size_t>(dof)];
                        if (std::isnan(value)) continue;

                        Vec6 constraint = Vec6::Constant(NAN);
                        constraint[dof] = value;

                        auto support = std::make_shared<bc::Support>(
                            support_region, constraint, support_orientation
                        );
                        if (*collector) {
                            (*collector)->add(std::move(support));
                        } else {
                            parser.select_collector_condition(bc::SUPPORT, std::move(support));
                        }
                    }
                }
            }
        );
    });
}

} // namespace fem::io::reader::commands