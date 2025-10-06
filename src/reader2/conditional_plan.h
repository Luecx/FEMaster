#pragma once
/**
 * @file conditional_plan.h
 * @brief Defines the ConditionalPlan class, representing a sequence of parsing segments
 *        that apply only when a specific condition is satisfied.
 *
 * @details
 * The ConditionalPlan acts as a fluent container for defining conditional parsing logic.
 * It links a Condition with one or more Segments and provides helper methods for
 * configuration and inspection. This is primarily used inside the command parsing
 * framework to define optional or context-sensitive blocks.
 *
 * @date 06.10.2025
 * @version 1.0
 * @author Finn Eggers
 */

#include <vector>
#include <string>
#include <initializer_list>
#include "condition.h"
#include "segment.h"

namespace fem::reader2 {

/**
 * @class ConditionalPlan
 * @brief Associates a Condition with a sequence of Segments.
 *
 * A ConditionalPlan defines a parsing plan that is only activated if its condition
 * evaluates to true. It supports fluent-style setup for chaining configuration calls.
 */
class ConditionalPlan {
private:
    Condition            _condition;     ///< Condition determining applicability
    std::string          _displayName;   ///< Human-readable plan name
    std::vector<Segment> _segments;      ///< Sequence of segments for this plan

public:
    /** @brief Default constructor initializing to an always-true condition. */
    ConditionalPlan();

    /** @brief Factory helper for fluent configuration. */
    static ConditionalPlan make();

    /** @brief Attach the condition that must hold for this plan. */
    ConditionalPlan& when(Condition condition);

    /** @brief Set a human-readable display name. */
    ConditionalPlan& name(std::string value);

    /** @brief Append a segment to the execution list. */
    ConditionalPlan& add(Segment segment);

    /** @brief Replace the segment list using an initializer list. */
    ConditionalPlan& segments(std::initializer_list<Segment> list);

    /** @brief Access the configured condition. */
    const Condition& cond() const;

    /** @brief Obtain the display name. */
    const std::string& display_name() const;

    /** @brief Inspect the configured segment sequence. */
    const std::vector<Segment>& segments() const;
};

} // namespace fem::reader2
