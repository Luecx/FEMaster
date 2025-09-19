/******************************************************************************
* @file timings.cpp
 * @brief Pretty-print compact timing tables for ConstraintBuilder stages.
 ******************************************************************************/

#include "timings.h"

#include "../../core/logging.h"
#include <iomanip>

namespace fem { namespace constraint {

void print_constraint_builder_timing(const std::vector<std::pair<std::string, Time>>& rows)
{
    logging::info(true, "");
    logging::info(true, "ConstraintBuilder timings (ms):");
    for (const auto& [name, t] : rows) {
        logging::info(true, "  ", std::left, std::setw(18), name, " ", std::right, std::setw(8), t);
    }
}

}} // namespace fem::constraint
