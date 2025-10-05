/**
 * @file constraint_report.cpp
 * @brief Implements utilities to format constraint equations for logging.
 *
 * Helper routines build readable strings from individual constraint equations
 * or full equation sets, including optional tagging of their origin.
 *
 * @see src/constraints/constraint_report.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_report.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace fem {
namespace constraint {
namespace {

/**
 * @brief Formats a single equation entry as `coeff * N<id>.<dof>`.
 */
std::string format_entry(const EquationEntry& entry, const EquationFormatOptions& opt) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(opt.precision);

    const Precision abs_coeff = std::abs(entry.coeff);
    if (abs_coeff != Precision(1) || abs_coeff == Precision(0)) {
        oss << abs_coeff << ' ';
    }
    oss << 'N' << entry.node_id << '.' << static_cast<int>(entry.dof);
    return oss.str();
}

/**
 * @brief Builds the right-hand side representation for an equation.
 */
std::string build_rhs(const Equation& equation, const EquationFormatOptions& opt) {
    std::ostringstream rhs;
    rhs.setf(std::ios::fixed, std::ios::floatfield);
    rhs.precision(opt.precision);
    const bool homogeneous = std::abs(equation.rhs) < Precision(1e-14);
    rhs << (homogeneous ? Precision(0) : equation.rhs);
    return rhs.str();
}

/**
 * @brief Builds the left-hand side representation for an equation.
 */
std::string build_lhs(const Equation& equation, const EquationFormatOptions& opt) {
    if (equation.entries.empty()) {
        return std::string("0");
    }

    std::ostringstream lhs;
    for (std::size_t idx = 0; idx < equation.entries.size(); ++idx) {
        const auto& entry = equation.entries[idx];
        const bool negative = std::signbit(entry.coeff);
        const bool first = idx == 0;

        if (first) {
            if (negative) {
                lhs << '-';
            }
        } else {
            lhs << (negative ? " - " : " + ");
        }
        lhs << format_entry(entry, opt);
    }
    return lhs.str();
}

/**
 * @brief Creates a compact tag describing the equation origin.
 */
std::string make_tag(EquationSourceKind source, Index idx) {
    char code = 0;
    switch (source) {
        case EquationSourceKind::Support:   code = 'S'; break;
        case EquationSourceKind::Connector: code = 'C'; break;
        case EquationSourceKind::Coupling:  code = 'U'; break;
        case EquationSourceKind::Tie:       code = 'T'; break;
        case EquationSourceKind::Manual:    code = 'M'; break;
        default:                            code = 0;  break;
    }
    if (code == 0) {
        return {};
    }
    return std::string(1, code) + std::to_string(idx);
}

} // namespace

/**
 * @copydoc format_equation
 */
std::string format_equation(const Equation& equation, const EquationFormatOptions& opt) {
    std::ostringstream line;
    const auto lhs = build_lhs(equation, opt);
    const auto rhs = build_rhs(equation, opt);
    line << lhs << " = " << rhs;
    return line.str();
}

/**
 * @copydoc format_equations
 */
std::vector<std::string> format_equations(const Equations& equations, const EquationFormatOptions& opt) {
    std::vector<std::string> lhs_parts;
    std::vector<std::string> rhs_parts;
    std::vector<std::string> tags;
    std::vector<bool> first_in_group;
    lhs_parts.reserve(equations.size());
    rhs_parts.reserve(equations.size());
    tags.reserve(equations.size());
    first_in_group.reserve(equations.size());

    std::size_t max_lhs = 0;
    std::size_t max_tag = 0;
    for (const auto& eq : equations) {
        auto lhs = build_lhs(eq, opt);
        auto rhs = build_rhs(eq, opt);
        max_lhs = std::max<std::size_t>(max_lhs, lhs.size());
        std::string tag = make_tag(eq.source, eq.source_index);
        if (!tag.empty()) {
            tag = '[' + tag + ']';
            max_tag = std::max<std::size_t>(max_tag, tag.size());
        }
        lhs_parts.emplace_back(std::move(lhs));
        rhs_parts.emplace_back(std::move(rhs));
        tags.emplace_back(std::move(tag));
    }

    std::vector<std::string> lines;
    lines.reserve(equations.size());
    for (std::size_t i = 0; i < tags.size(); ++i) {
        if (i == 0) {
            first_in_group.push_back(true);
        } else {
            first_in_group.push_back(tags[i] != tags[i - 1]);
        }
    }

    for (std::size_t i = 0; i < equations.size(); ++i) {
        std::ostringstream oss;
        oss.setf(std::ios::fixed, std::ios::floatfield);
        if (max_tag > 0) {
            if (first_in_group[i] || tags[i].empty()) {
                oss << std::left << std::setw(static_cast<int>(max_tag)) << tags[i] << ' ';
            } else {
                oss << std::string(static_cast<std::size_t>(max_tag), ' ') << ' ';
            }
        }
        oss << std::left << std::setw(static_cast<int>(max_lhs)) << lhs_parts[i]
            << " = " << rhs_parts[i];
        lines.emplace_back(oss.str());
    }
    return lines;
}

} // namespace constraint
} // namespace fem
