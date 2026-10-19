/**
 * @file condition_state.h
 * @brief Defines step-owned structural and thermal condition sets.
 *
 * A ConditionState groups physical condition definitions by input family.
 * The parser resolves input identifiers and OP=MOD/NEW against these groups;
 * load cases own independent initial and target states. The finite-element
 * model owns neither active prescriptions nor step-history membership.
 *
 * @see Condition
 * @see LoadCase
 */

#pragma once

#include "condition.h"

#include <array>
#include <cstddef>
#include <unordered_set>
#include <utility>

namespace fem::bc {

/**
 * @brief Independent input-history scope for a physical condition.
 *
 * Families retain the distinction between keyword operations without
 * imposing a physical interpolation rule. Load evaluation distinguishes
 * structural loads from support and thermal prescriptions separately.
 */
enum ConditionFamily : std::size_t {
    SUPPORT = 0,
    CLOAD,
    DLOAD,
    DSLOAD,
    PLOAD,
    VLOAD,
    INERTIAL_LOAD,
    TEMPERATURE,
    HEAT_FLUX,
    CONVECTION,
    N_CONDITION_FAMILIES
};

/**
 * @brief Grouped prescriptions representing one step boundary state.
 *
 * A load case holds a frozen initial state and an independently editable
 * target state. Copying the container shares condition definitions until a
 * replacement is introduced, whereas clone() creates independent physical
 * objects using their polymorphic copy operations. Shared regions and named
 * amplitudes remain reusable immutable definitions.
 *
 * Pointer identity prevents inserting a single definition twice into one
 * family; distinct overlapping prescriptions remain additive.
 */
struct ConditionState {
    using Conditions = std::unordered_set<Condition::Ptr>;

    std::array<Conditions, N_CONDITION_FAMILIES> conditions;

    // Membership and access grouped by input-history family
    bool add(ConditionFamily family, Condition::Ptr condition) {
        return conditions[family].insert(std::move(condition)).second;
    }

    bool remove(ConditionFamily family, const Condition::Ptr& condition) {
        return conditions[family].erase(condition) != 0;
    }

    bool contains(ConditionFamily family, const Condition::Ptr& condition) const {
        return conditions[family].find(condition) != conditions[family].end();
    }

    const Conditions& get(ConditionFamily family) const {
        return conditions[family];
    }

    // Deep-copy physical condition objects while retaining shared immutable definitions
    ConditionState clone() const {
        ConditionState result;
        for (std::size_t i = 0; i < N_CONDITION_FAMILIES; ++i) {
            for (const auto& condition : conditions[i]) {
                result.conditions[i].insert(condition->clone());
            }
        }
        return result;
    }
};

} // namespace fem::bc
