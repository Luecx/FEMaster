#pragma once
/**
 * @file line_range.h
 * @brief Declares Range, inclusive bounds limiting how many data groups a segment may consume.
 */

#include <cstddef>
#include <limits>

namespace fem::reader2 {

    /**
     * @brief Inclusive bounds that limit how many data groups a segment may consume.
     */
    struct Range {
        size_t _min = 0;
        size_t _max = std::numeric_limits<size_t>::max();

        /** @brief Factory for fluent-style construction. */
        static Range make();

        /** @brief Set the lower bound. */
        Range& min(size_t value);

        /** @brief Set the upper bound. */
        Range& max(size_t value);

        /** @brief Retrieve the lower bound. */
        size_t min() const;

        /** @brief Retrieve the upper bound. */
        size_t max() const;
    };

} // namespace fem::reader2
