/******************************************************************************
 * @file timings.h
 * @brief Declares helpers to print builder timing summaries.
 *
 * Utilities in this module produce concise tables summarising the durations of
 * each constraint-building stage.
 *
 * @see src/constraints/builder/timings.cpp
 * @see src/constraints/builder/builder.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../../core/timer.h"

#include <string>
#include <utility>
#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @brief Prints a formatted table of timing measurements.
 *
 * @param rows Name/time pairs that describe each pipeline stage.
 ******************************************************************************/
void print_constraint_builder_timing(const std::vector<std::pair<std::string, Time>>& rows);

} // namespace constraint
} // namespace fem
