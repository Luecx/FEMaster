/**
 * @file constraint_report.h
 * @brief Declares helpers to format constraint equations for logging.
 *
 * Utility functions are provided to turn individual equations or full sets
 * into readable text lines that aid debugging and reporting.
 *
 * @see src/constraints/constraint_report.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "equation.h"

#include <string>
#include <vector>

namespace fem {
namespace constraint {

/**
 * @struct EquationFormatOptions
 * @brief Configures the formatting of constraint equations.
 */
struct EquationFormatOptions {
    int precision = 3;   ///< Decimal precision used for coefficients.
    bool show_sign = true; ///< Whether to prefix positive coefficients with a plus sign.
};

/**
 * @brief Formats a single constraint equation into a human-readable string.
 *
 * @param equation Equation to format.
 * @param opt Formatting options.
 * @return std::string Generated textual representation.
 */
std::string format_equation(const Equation& equation, const EquationFormatOptions& opt = {});

/**
 * @brief Formats an entire collection of equations line-by-line.
 *
 * @param equations Equation list to format.
 * @param opt Formatting options.
 * @return std::vector<std::string> Textual representations for all equations.
 */
std::vector<std::string> format_equations(const Equations& equations, const EquationFormatOptions& opt = {});

} // namespace constraint
} // namespace fem
