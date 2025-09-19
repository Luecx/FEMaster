/******************************************************************************
* @file timings.h
 * @brief Pretty-print compact timing tables for ConstraintBuilder stages.
 ******************************************************************************/

#pragma once
#include "../../core/timer.h"
#include <string>
#include <vector>
#include <utility>

namespace fem { namespace constraint {

void print_constraint_builder_timing(const std::vector<std::pair<std::string, Time>>& rows);

}} // namespace
