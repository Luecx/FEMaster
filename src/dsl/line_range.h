/**
* @file line_range.h
 * @brief Declares a small helper that constrains how many data lines a segment may consume.
 *
 * `LineRange` is used by the DSL to describe the **inclusive** minimum/maximum number of
 * consecutive data lines that a segment expects. It is purely a configuration object; the
 * engine enforces these bounds during parsing.
 *
 * - Default is `[1, +inf)` meaning: at least one data line, unbounded maximum.
 * - Both bounds are **inclusive**.
 * - Use the fluent setters `min(...)` and `max(...)` for readability:
 *   `Segment::make().range(LineRange{}.min(1).max(3))`
 *
 * @see segment.h
 * @see pattern.h
 *
 * @author
 *   Finn Eggers (project)
 * @date
 *   12.10.2025
 */

#pragma once
#include <cstddef>
#include <limits>

namespace fem {
namespace dsl {

/**
 * @struct LineRange
 * @brief Inclusive lower/upper bounds on the number of data lines a segment may read.
 */
struct LineRange {
    /** Inclusive lower bound (defaults to 1). */
    std::size_t min_ = 1;

    /** Inclusive upper bound (defaults to the maximum representable `size_t`). */
    std::size_t max_ = std::numeric_limits<std::size_t>::max();

    /**
     * @brief Sets the inclusive minimum number of lines.
     * @param v New minimum (must be ≤ current `max_` to be meaningful).
     * @return Reference to `*this` for fluent chaining.
     */
    LineRange& min(std::size_t v) {
        min_ = v;
        return *this;
    }

    /**
     * @brief Sets the inclusive maximum number of lines.
     * @param v New maximum (should be ≥ current `min_` to be meaningful).
     * @return Reference to `*this` for fluent chaining.
     */
    LineRange& max(std::size_t v) {
        max_ = v;
        return *this;
    }
};

} // namespace dsl
} // namespace fem
