/**
 * @file variant.h
 * @brief Declares `Variant`, a selectable layout of segments under a key/parent condition.
 *
 * A `Variant` groups one or more `Segment`s and is optionally guarded by a key/parent
 * condition. During parsing, the engine selects the **first** variant whose condition
 * evaluates to true (or the first without a condition).
 *
 * Example:
 * @code
 * Variant::make()
 *   .when( all_of({ key_equals("MODE", {"ORTHO"}), parent_is("MATERIAL") }) )
 *   .segment( consecutive data block #1 )
 *   .segment( consecutive data block #2 );
 * @endcode
 *
 * Notes:
 *  - `when(...)` stores the provided condition by value.
 *  - This implementation enforces a simple "no-overlap" style:
 *      * All previously added segments of the same variant must be single-line
 *        (`pattern.is_multiline() == false`) and have a fixed line count
 *        (`range.min_ == range.max_`). This guarantees that each segment consumes
 *        a deterministic, non-extensible chunk, so following segments cannot be
 *        "überlesen".
 *      * Additionally, each segment must require at least one token.
 *    If you need multi-line data, place it in the **last** segment of a variant,
 *    or split layouts into separate variants with `when(...)` guards.
 *
 * @see condition.h
 * @see segment.h
 * @see command.h
 *
 * @date 12.10.2025
 * @author
 *   Finn Eggers (project)
 */

#pragma once
#include <utility>
#include <vector>
#include <stdexcept>
#include <sstream>
#include "segment.h"
#include "condition.h"

namespace fem {
namespace dsl {

/**
 * @class Variant
 * @brief A conditional collection of segments representing one admissible data layout.
 *
 * A command may define multiple variants; the engine will pick the first variant whose
 * condition is satisfied (or has no condition). Each variant contains an ordered list of
 * segments that are executed sequentially.
 */
struct Variant {
    /** Whether a condition has been explicitly provided via `when(...)`. */
    bool has_condition_ = false;

    /** Stored condition; meaningful only if `has_condition_` is `true`. */
    Condition condition_{};

    /** Ordered list of segments that compose this variant. */
    std::vector<Segment> _segments;

    /** Optional short description used by documentation printers. */
    std::string doc_;

    /**
     * @brief Factory for fluent construction.
     */
    static Variant make() { return {}; }

    /**
     * @brief Sets the condition that must hold for this variant to be chosen.
     *
     * @param c Condition tree to store by value.
     * @return Reference to `*this` for fluent chaining.
     */
    Variant& when(Condition c) {
        condition_ = std::move(c);
        has_condition_ = true;
        return *this;
    }

    /**
     * @brief Appends a segment to the variant.
     *
     * Enforces a conservative "no-overlap" policy for already added segments:
     *  - All previously added segments must be single-line (`!pattern.is_multiline()`).
     *  - All previously added segments must have a fixed line count (`range.min_ == range.max_`).
     *  - The new segment's pattern must require at least one token.
     *
     * Rationale: Multi-line segments or unbounded ranges in the middle of a variant
     * could consume an unpredictable number of lines/tokens, making nachfolgende
     * Segmente fragil. Wenn du Multi-line brauchst, platziere es als **letztes**
     * Segment oder verwende getrennte Variants.
     *
     * @param s Segment to append (moved in).
     * @throws std::runtime_error if style rules are violated.
     * @return Reference to `*this` for fluent chaining.
     */
    Variant& segment(Segment s) {
        // Guard: the new segment itself must require tokens (sanity).
        if (s._pattern.required_tokens() == 0) {
            throw std::runtime_error("Variant::segment(): pattern requires 0 tokens (not allowed).");
        }

        // Check all previously added segments for "no-overlap" constraints.
        for (std::size_t i = 0; i < _segments.size(); ++i) {
            const Segment& prev = _segments[i];

            // previous must not be multiline
            if (prev._pattern.is_multiline()) {
                std::ostringstream os;
                os << "Variant::segment(): previous segment #" << (i+1)
                   << " is multiline; only the last segment of a variant may be multiline.";
                throw std::runtime_error(os.str());
            }

            // previous must have fixed number of lines (min == max)
            if (prev._range.min_ != prev._range.max_) {
                std::ostringstream os;
                os << "Variant::segment(): previous segment #" << (i+1)
                   << " has non-fixed line range [" << prev._range.min_ << ".." << prev._range.max_
                   << "]; only fixed (min==max) allowed before the last segment.";
                throw std::runtime_error(os.str());
            }

            // previous must also require at least one token
            if (prev._pattern.required_tokens() == 0) {
                std::ostringstream os;
                os << "Variant::segment(): previous segment #" << (i+1)
                   << " has a pattern requiring 0 tokens.";
                throw std::runtime_error(os.str());
            }
        }

        _segments.push_back(std::move(s));
        return *this;
    }

    /**
     * @brief Sets a human-readable description for this variant (shown in help output).
     * @param d Short description text.
     * @return Reference to `*this` for fluent chaining.
     */
    Variant& doc(std::string d) {
        doc_ = std::move(d);
        return *this;
    }
};

} // namespace dsl
} // namespace fem
