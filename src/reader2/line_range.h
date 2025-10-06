#pragma once
#include <cstddef>
#include <limits>

namespace fem::reader2 {

/**
 * \brief Inclusive bounds that limit how many data groups a segment may consume.
 */
struct Range {
    size_t _min = 0;
    size_t _max = std::numeric_limits<size_t>::max();

    /**
     * \brief Factory for fluent-style construction.
     */
    static Range make()
    {
        return Range{};
    }

    /**
     * \brief Set the lower bound.
     */
    Range& min(size_t value)
    {
        _min = value;
        return *this;
    }

    /**
     * \brief Set the upper bound.
     */
    Range& max(size_t value)
    {
        _max = value;
        return *this;
    }

    /**
     * \brief Retrieve the lower bound.
     */
    size_t min() const
    {
        return _min;
    }

    /**
     * \brief Retrieve the upper bound.
     */
    size_t max() const
    {
        return _max;
    }
};

} // namespace fem::reader2
