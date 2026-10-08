/**
 * @file register_loadcase_newmark.cpp
 * @brief Registers Newmark-beta integration parameters for transient analyses.
 *
 * The `NEWMARK` child command reads beta and gamma for the implicit fixed-step
 * time integrator and stores them on the active `Transient` load case. Default
 * values remain beta = 0.25 and gamma = 0.5 when the command is absent.
 *
 * Time stepping, effective operator construction and state advancement remain
 * responsibilities of the transient solver.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>

#include "../parser.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"

#include "../../../loadcase/linear_transient.h"

namespace fem::io::reader::commands {

void register_loadcase_newmark(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("NEWMARK", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Set Newmark-β integration parameters (β, γ). Defaults are 0.25, 0.5.");

        command.variant(
            fem::io::dsl::Variant::make()
                .doc("One data line: β, γ")
                .data(
                    fem::io::dsl::Pattern::make()
                                .fixed<fem::Precision, 2>("NEWMARK", "β, γ parameters for Newmark-β.")
                                .defaults(fem::Precision{0}),
                    [&parser](const std::array<fem::Precision, 2>& bg) {
                            auto* base = parser.active_loadcase();
                            logging::error(base != nullptr, "NEWMARK must appear inside *LOADCASE.");

                            if (auto* lc = base->as<fem::loadcase::Transient>()) {
                                lc->set_newmark(bg[0], bg[1]);
                                return;
                            }

                            logging::error(false, "NEWMARK not supported for loadcase type " + base->type_name());
                        },
                    fem::io::dsl::LineRange{}.min(1).max(1)
                )
        );
    });
}

} // namespace fem::io::reader::commands
