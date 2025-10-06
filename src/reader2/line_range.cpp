/**
* @file line_range.cpp
 * @brief Implements Range, inclusive bounds for segment group counts.
 */

#include "line_range.h"

namespace fem::reader2 {

Range Range::make() {
    return Range{};
}

Range& Range::min(size_t value) {
    _min = value;
    return *this;
}

Range& Range::max(size_t value) {
    _max = value;
    return *this;
}

size_t Range::min() const {
    return _min;
}

size_t Range::max() const {
    return _max;
}

} // namespace fem::reader2
