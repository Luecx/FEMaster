/**
 * @file register_specific_heat.cpp
 * @brief Reads constant specific heat for material definitions.
 *
 * This command stores a temperature-independent specific heat on the active
 * material. The value is available for future transient thermal analyses;
 * registering it does not change any current analysis procedure.
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"

#include <cmath>

namespace fem::io::reader::commands {

/**
 * Registers the Abaqus `*SPECIFIC HEAT` material keyword.
 *
 * Only one constant scalar is accepted. Temperature-dependent tables are not
 * interpreted as constant data and remain unsupported.
 */
void register_specific_heat(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SPECIFICHEAT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign constant specific heat to the active material.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<Precision>().name("CP").desc("Constant specific heat")
                )
                .bind([&model](Precision specific_heat) {
                    logging::error(std::isfinite(specific_heat) && specific_heat > Precision(0),
                        "SPECIFIC HEAT requires a finite, positive coefficient");
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr,
                        "SPECIFIC HEAT requires an active material context");
                    material->set_thermal_specific_heat(specific_heat);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
