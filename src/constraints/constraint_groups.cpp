/**
 * @file constraint_groups.cpp
 * @brief Implements logging helpers for grouped constraint equations.
 *
 * @see src/constraints/constraint_groups.h
 */

#include "constraint_groups.h"

#include "constraint_report.h"
#include "../core/logging.h"

namespace fem {
namespace constraint {

void ConstraintGroups::report() const {
    using fem::logging::info;

    const constraint::EquationFormatOptions format{};
    const struct {
        const char* name;
        const constraint::Equations& equations;
    } categories[] = {
        {"Supports", supports},
        {"Connectors", connectors},
        {"Couplings", couplings},
        {"Ties", ties},
        {"RBMs", rbms},
        {"Other", others},
    };

    info(true, "Constraint overview:");
    logging::up();
    for (const auto& category : categories) {
        info(true, category.name, " (", category.equations.size(), ")");
        if (!category.equations.empty()) {
            logging::up();
            for (const auto& line : constraint::format_equations(category.equations, format)) {
                info(true, "    ", line);
            }
            logging::down();
        }
    }
    logging::down();
}

} // namespace constraint
} // namespace fem
