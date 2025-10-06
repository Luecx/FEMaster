#pragma once
#include <string>
#include "line_range.h"
#include "pattern.h"

namespace fem::reader2 {

/**
 * \brief Associates a Pattern with multiplicity and presentation metadata.
 */
struct Segment {
private:
    Range       _range;
    Pattern     _pattern;
    std::string _label;

public:
    static Segment make()
    {
        return Segment{};
    }

    Segment& range(Range value)
    {
        _range = std::move(value);
        return *this;
    }

    Segment& pattern(Pattern value)
    {
        _pattern = std::move(value);
        return *this;
    }

    Segment& label(std::string value)
    {
        _label = std::move(value);
        return *this;
    }

    const Range& range() const
    {
        return _range;
    }

    const Pattern& pattern() const
    {
        return _pattern;
    }

    const std::string& label() const
    {
        return _label;
    }
};

} // namespace fem::reader2
