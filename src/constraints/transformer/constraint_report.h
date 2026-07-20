/**
 * @file constraint_report.h
 * @brief Declares textual formatting utilities for linear constraint equations.
 *
 * The functions declared in this file convert individual constraint equations
 * and complete equation collections into human-readable strings. The generated
 * output is intended for diagnostic logging, debugging and inspection of the
 * assembled constraint definitions before they are transformed into the global
 * sparse system
 *
 *     C u = d.
 *
 * Each equation is represented using nodal identifiers and local degree-of-
 * freedom indices. Collections can additionally be annotated with compact tags
 * describing the source of each equation and aligned into a readable tabular
 * layout.
 *
 * @see constraint_report.cpp
 * @see ../types/equation.h
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../types/equation.h"

#include <string>
#include <vector>

namespace fem {
namespace constraint {

/**
 * @brief Configures the textual representation of constraint equations.
 *
 * The options control the numerical formatting used when coefficients and
 * right-hand-side values are converted to text. They do not affect the
 * underlying equation data.
 */
struct EquationFormatOptions {
    // Number of digits written after the decimal point for coefficients and
    // right-hand-side values. Formatting uses fixed-point notation.
    int precision = 3;

    // Indicates whether positive terms should be shown with an explicit leading
    // plus sign. This option is part of the public formatting interface, but the
    // current implementation determines signs from the term position and does
    // not evaluate this value directly.
    bool show_sign = true;
};

// Convert one linear constraint equation into a compact readable expression of
// the form
//
//     coefficient * N<node>.<dof> + ... = rhs.
//
// Unit coefficients are omitted, negative signs are emitted separately and an
// empty left-hand side is represented by zero.
std::string format_equation(const Equation& equation, const EquationFormatOptions& opt = {});

// Convert a complete equation collection into aligned text lines. Equations are
// annotated with compact source tags where available, consecutive equations
// from the same source group suppress repeated tags and all left-hand sides are
// padded to a common width.
std::vector<std::string> format_equations(const Equations& equations, const EquationFormatOptions& opt = {});

} // namespace constraint
} // namespace fem
