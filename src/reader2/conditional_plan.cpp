/**
* @file conditional_plan.cpp
 * @brief Implements the ConditionalPlan class for conditional parsing sequences.
 *
 * @date 06.10.2025
 * @version 1.0
 * @see conditional_plan.h
 */

#include "conditional_plan.h"

namespace fem::reader2 {

    ConditionalPlan::ConditionalPlan()
        : _condition(Condition::always()), _displayName(), _segments() {}

    ConditionalPlan ConditionalPlan::make() {
        return ConditionalPlan{};
    }

    ConditionalPlan& ConditionalPlan::when(Condition condition) {
        _condition = std::move(condition);
        return *this;
    }

    ConditionalPlan& ConditionalPlan::name(std::string value) {
        _displayName = std::move(value);
        return *this;
    }

    ConditionalPlan& ConditionalPlan::add(Segment segment) {
        _segments.push_back(std::move(segment));
        return *this;
    }

    ConditionalPlan& ConditionalPlan::segments(std::initializer_list<Segment> list) {
        _segments.assign(list.begin(), list.end());
        return *this;
    }

    const Condition& ConditionalPlan::cond() const {
        return _condition;
    }

    const std::string& ConditionalPlan::display_name() const {
        return _displayName;
    }

    const std::vector<Segment>& ConditionalPlan::segments() const {
        return _segments;
    }

} // namespace fem::reader2
