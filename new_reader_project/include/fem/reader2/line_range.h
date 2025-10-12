
#pragma once
#include <cstddef>
#include <limits>

namespace fem { namespace reader2 {
struct LineRange {
    std::size_t min_ = 1;
    std::size_t max_ = std::numeric_limits<std::size_t>::max();
    LineRange& min(std::size_t v){ min_ = v; return *this; }
    LineRange& max(std::size_t v){ max_ = v; return *this; }
};
}} // ns
