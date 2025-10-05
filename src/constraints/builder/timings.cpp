/**
 * @file timings.cpp
 * @brief Prints compact timing tables for constraint-building stages.
 *
 * The helper formats timing data collected during the builder pipeline for
 * easier inspection.
 *
 * @see src/constraints/builder/timings.h
 * @see src/constraints/builder/builder.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "timings.h"

#include "../../core/logging.h"

#include <iomanip>

namespace fem {
namespace constraint {

/**
 * @copydoc print_constraint_builder_timing
 */
void print_constraint_builder_timing(const std::vector<std::pair<std::string, Time>>& rows) {
    logging::info(true, "");
    logging::info(true, "ConstraintBuilder timings (ms):");
    for (const auto& [name, duration] : rows) {
        logging::info(true,
                      "  ", std::left, std::setw(18), name,
                      " ", std::right, std::setw(8), duration);
    }
}

} // namespace constraint
} // namespace fem
