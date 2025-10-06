#pragma once
#include <vector>
#include <string>
#include <initializer_list>
#include "condition.h"
#include "segment.h"

namespace fem::reader2 {

/**
 * \brief Describes a sequence of segments that apply when a condition is satisfied.
 */
struct ConditionalPlan {
private:
    Condition            _condition = Condition::always();
    std::string          _displayName;
    std::vector<Segment> _segments;

public:
    /**
     * \brief Factory helper for fluent configuration.
     */
    static ConditionalPlan make()
    {
        return ConditionalPlan{};
    }

    /**
     * \brief Attach the condition that needs to hold for this plan.
     */
    ConditionalPlan& when(Condition condition)
    {
        _condition = std::move(condition);
        return *this;
    }

    /**
     * \brief Set a human-friendly display name.
     */
    ConditionalPlan& name(std::string value)
    {
        _displayName = std::move(value);
        return *this;
    }

    /**
     * \brief Append a segment to the execution list.
     */
    ConditionalPlan& add(Segment segment)
    {
        _segments.push_back(std::move(segment));
        return *this;
    }

    /**
     * \brief Replace the segment list using an initializer list.
     */
    ConditionalPlan& segments(std::initializer_list<Segment> list)
    {
        _segments.assign(list.begin(), list.end());
        return *this;
    }

    /**
     * \brief Access the configured condition.
     */
    const Condition& cond() const
    {
        return _condition;
    }

    /**
     * \brief Obtain the display name.
     */
    const std::string& display_name() const
    {
        return _displayName;
    }

    /**
     * \brief Inspect the configured segment sequence.
     */
    const std::vector<Segment>& segments() const
    {
        return _segments;
    }
};

} // namespace fem::reader2
