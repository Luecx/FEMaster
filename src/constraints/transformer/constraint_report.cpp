
/**
 * @file constraint_report.cpp
 * @brief Implements textual formatting of linear constraint equations.
 *
 * The formatting utilities convert equations of the form
 *
 *     sum_i c_i u_i = d
 *
 * into readable expressions using nodal identifiers and degree-of-freedom
 * indices. Individual equations can be formatted independently, while complete
 * equation collections are aligned and optionally annotated with compact tags
 * that identify the origin of each equation.
 *
 * The implementation is deliberately independent of sparse constraint-system
 * assembly. It operates directly on the symbolic `Equation` representation and
 * is therefore suitable for diagnostics before inactive DOFs are removed or
 * equations are transformed into the global matrix `C`.
 *
 * @see constraint_report.h
 * @see ../types/equation.h
 *
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

// Format one term of a constraint equation as
//
//     <absolute coefficient> N<node id>.<dof>
//
// The algebraic sign is intentionally not written here. It is handled by
// build_lhs(), which can distinguish the first term from subsequent terms and
// insert either a leading minus sign or an infix plus/minus separator.
std::string format_entry(const EquationEntry& entry, const EquationFormatOptions& opt) {
    // Use an isolated string stream so the numerical formatting of one
    // coefficient does not affect the surrounding equation stream.
    std::ostringstream oss;

    // Write coefficients in fixed-point notation using the user-selected
    // number of decimal places.
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(opt.precision);

    // Signs are emitted by build_lhs(), so only the coefficient magnitude is
    // represented in the term itself.
    const Precision abs_coeff = std::abs(entry.coeff);

    // Omit coefficients with unit magnitude to produce expressions such as
    //
    //     N2.0 - N5.0
    //
    // instead of
    //
    //     1.000 N2.0 - 1.000 N5.0.
    //
    // The explicit zero check ensures that a zero coefficient remains visible
    // if such an entry is present in the symbolic equation.
    if (abs_coeff != Precision(1) || abs_coeff == Precision(0)) {
        oss << abs_coeff << ' ';
    }

    // Encode the addressed degree of freedom through its node identifier and
    // local DOF index.
    oss << 'N' << entry.node_id << '.' << static_cast<int>(entry.dof);

    return oss.str();
}

// Format the scalar right-hand side of an equation in fixed-point notation.
std::string build_rhs(const Equation& equation, const EquationFormatOptions& opt) {
    std::ostringstream rhs;

    // Apply the same numerical representation used for equation coefficients.
    rhs.setf(std::ios::fixed, std::ios::floatfield);
    rhs.precision(opt.precision);

    // Treat very small right-hand-side values as exactly zero in the textual
    // representation. This suppresses numerical noise and makes homogeneous
    // equations immediately recognizable in diagnostic output.
    const bool homogeneous = std::abs(equation.rhs) < Precision(1e-14);

    // Preserve the original value for genuinely inhomogeneous equations and
    // emit an exact zero otherwise.
    rhs << (homogeneous ? Precision(0) : equation.rhs);

    return rhs.str();
}

// Build the complete left-hand side of one equation while preserving the sign
// and order of every symbolic entry.
std::string build_lhs(const Equation& equation, const EquationFormatOptions& opt) {
    // An equation without entries has a mathematically empty sum and is shown
    // explicitly as zero.
    if (equation.entries.empty()) {
        return std::string("0");
    }

    std::ostringstream lhs;

    // Process entries in their stored order. No sorting or coefficient
    // consolidation is performed by the reporting layer.
    for (std::size_t idx = 0; idx < equation.entries.size(); ++idx) {
        const auto& entry = equation.entries[idx];

        // signbit() also detects a negative signed zero, preserving the sign
        // carried by the stored floating-point coefficient.
        const bool negative = std::signbit(entry.coeff);
        const bool first = idx == 0;

        // The first term receives no leading plus sign. A negative first term
        // is prefixed only by a minus sign.
        if (first) {
            if (negative) {
                lhs << '-';
            }
        } else {
            // Subsequent terms use an infix separator that contains both the
            // algebraic sign and surrounding whitespace.
            lhs << (negative ? " - " : " + ");
        }

        // Append the unsigned term representation. The sign has already been
        // written by the logic above.
        lhs << format_entry(entry, opt);
    }

    return lhs.str();
}

// Create the compact source identifier used to group equations in multi-line
// diagnostic output.
//
// The generated tag consists of one source-specific letter followed by the
// source index, for example `S3` for a support equation or `T1` for a tie
// equation.
std::string make_tag(EquationSourceKind source, Index idx) {
    char code = 0;

    // Map each known equation category to a short, visually distinct prefix.
    switch (source) {
        case EquationSourceKind::Support:   code = 'S'; break;
        case EquationSourceKind::Connector: code = 'C'; break;
        case EquationSourceKind::Coupling:  code = 'U'; break;
        case EquationSourceKind::Tie:       code = 'T'; break;
        case EquationSourceKind::Rbm:       code = 'R'; break;
        case EquationSourceKind::Manual:    code = 'M'; break;
        default:                            code = 0;  break;
    }

    // Unknown or intentionally untagged equation sources produce no visible
    // source identifier.
    if (code == 0) {
        return {};
    }

    // Combine the source prefix with its source-local index.
    return std::string(1, code) + std::to_string(idx);
}

} // namespace

// Format one equation as a complete single-line expression
//
//     <left-hand side> = <right-hand side>.
std::string format_equation(const Equation& equation, const EquationFormatOptions& opt) {
    std::ostringstream line;

    // Build both sides independently so their specialized formatting remains
    // encapsulated in the helper functions above.
    const auto lhs = build_lhs(equation, opt);
    const auto rhs = build_rhs(equation, opt);

    // Join both representations using the conventional equality separator.
    line << lhs << " = " << rhs;

    return line.str();
}

// Format an equation collection as aligned diagnostic lines with optional
// source tags and visual grouping of consecutive equations from the same
// source.
std::vector<std::string> format_equations(const Equations& equations, const EquationFormatOptions& opt) {
    // Store the separately formatted components before producing final lines.
    // This first pass is necessary to determine the maximum widths used for
    // column alignment.
    std::vector<std::string> lhs_parts;
    std::vector<std::string> rhs_parts;
    std::vector<std::string> tags;
    std::vector<bool> first_in_group;

    // Reserve one entry per equation to avoid repeated reallocations while the
    // formatted components are collected.
    lhs_parts.reserve(equations.size());
    rhs_parts.reserve(equations.size());
    tags.reserve(equations.size());
    first_in_group.reserve(equations.size());

    // Track the widest left-hand side and source tag. These widths define the
    // alignment columns of the final report.
    std::size_t max_lhs = 0;
    std::size_t max_tag = 0;

    // Pre-format every equation and collect the maximum field widths.
    for (const auto& eq : equations) {
        auto lhs = build_lhs(eq, opt);
        auto rhs = build_rhs(eq, opt);

        // Expand the left-hand-side column width whenever a longer expression
        // is encountered.
        max_lhs = std::max<std::size_t>(max_lhs, lhs.size());

        // Convert the equation origin into a compact source tag.
        std::string tag = make_tag(eq.source, eq.source_index);

        if (!tag.empty()) {
            // Enclose visible source tags in brackets and update the required
            // tag-column width.
            tag = '[' + tag + ']';
            max_tag = std::max<std::size_t>(max_tag, tag.size());
        }

        // Move the completed strings into their storage arrays because the
        // local temporaries are no longer needed.
        lhs_parts.emplace_back(std::move(lhs));
        rhs_parts.emplace_back(std::move(rhs));
        tags.emplace_back(std::move(tag));
    }

    std::vector<std::string> lines;
    lines.reserve(equations.size());

    // Mark the first equation of every consecutive source group. Repeated tags
    // inside one group are suppressed later to make related equations easier to
    // scan visually.
    for (std::size_t i = 0; i < tags.size(); ++i) {
        if (i == 0) {
            // The first equation always starts a new group.
            first_in_group.push_back(true);
        } else {
            // A group boundary occurs whenever the current tag differs from the
            // immediately preceding one.
            first_in_group.push_back(tags[i] != tags[i - 1]);
        }
    }

    // Assemble the final aligned output one equation at a time.
    for (std::size_t i = 0; i < equations.size(); ++i) {
        std::ostringstream oss;

        // Keep fixed-point formatting active for consistency, although the
        // numerical components have already been converted to strings.
        oss.setf(std::ios::fixed, std::ios::floatfield);

        // Emit the optional source-tag column only when at least one equation
        // in the collection has a visible tag.
        if (max_tag > 0) {
            if (first_in_group[i] || tags[i].empty()) {
                // Show the tag at the beginning of a new group. Empty tags are
                // also written explicitly as padding so the equation columns
                // remain aligned.
                oss << std::left << std::setw(static_cast<int>(max_tag)) << tags[i] << ' ';
            } else {
                // Suppress repeated tags within a consecutive group while
                // preserving exactly the same horizontal space.
                oss << std::string(static_cast<std::size_t>(max_tag), ' ') << ' ';
            }
        }

        // Left-align every left-hand side to the maximum detected width, then
        // append the equality sign and the equation-specific right-hand side.
        oss << std::left << std::setw(static_cast<int>(max_lhs)) << lhs_parts[i]
            << " = " << rhs_parts[i];

        lines.emplace_back(oss.str());
    }

    return lines;
}

} // namespace constraint
} // namespace fem

