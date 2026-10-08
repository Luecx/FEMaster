/**
 * @file condition_manager.cpp
 * @brief Implements identity-based ownership of active model conditions.
 *
 * Each operation acts on one input-history family of the model-owned
 * ConditionManager. A condition's membership is determined solely by its
 * shared pointer identity. Neither the physical target of a condition nor
 * its magnitudes and amplitudes are inspected in this implementation.
 *
 * The input reader resolves keyword-specific history, OP=MOD/OP=NEW and
 * original input identifiers. Model assembly subsequently reads these sets
 * and invokes the physical Condition::apply() implementations. Keeping
 * these responsibilities separate avoids coupling the finite-element
 * assembly layer to any particular input format.
 *
 * @see ConditionManager
 * @see Condition
 *
 * @author Finn Eggers
 * @date 08.10.2026
 */

#include "condition_manager.h"

#include "../core/logging.h"

#include <utility>

namespace fem::bc {

/**
 * @brief Adds one condition reference to an active history family.
 *
 * Hashing and equality operate on the stored shared pointer, so the operation
 * cannot collapse two distinct conditions merely because they prescribe the
 * same region or have equal numerical values. The same object may also appear
 * in different families because each set is independent.
 *
 * This method does not interpret input identifiers or replace semantic targets.
 * The input reader resolves those operations before calling add().
 *
 * @param family History family receiving the reference.
 * @param condition Non-null physical definition to insert.
 * @return True only if the pointer was not already present in this family.
 */
bool ConditionManager::add(ConditionFamily family, Condition::Ptr condition) {
    // Refuse invalid definitions before adding them to active model storage.
    logging::error(condition != nullptr,
        "ConditionManager cannot add a null condition");

    // Preserve distinct physical objects while suppressing duplicate pointers.
    return conditions_[family].insert(std::move(condition)).second;
}

/**
 * @brief Removes one exact condition pointer from an active family.
 *
 * The operation does not search for an equivalent input target. Removing a
 * manager reference does not destroy the definition while the parser, a named
 * collector or another owner still holds a shared pointer to it.
 *
 * @param family Family from which the reference is removed.
 * @param condition Pointer identity to remove.
 * @return True if the pointer was present and has been erased.
 */
bool ConditionManager::remove(ConditionFamily family, const Condition::Ptr& condition) {
    return conditions_[family].erase(condition) != 0;
}

/**
 * @brief Tests whether one exact condition object is active in a family.
 *
 * This is an identity and membership query, not a semantic comparison. In
 * particular, two distinct objects that address the same node or surface
 * are not interchangeable, even when Condition::matches() would accept them.
 *
 * @param family Family being queried.
 * @param condition Pointer identity being queried.
 * @return True only when this exact pointer occurs in the selected family.
 */
bool ConditionManager::contains(ConditionFamily family, const Condition::Ptr& condition) const {
    return conditions_[family].find(condition) != conditions_[family].end();
}

/**
 * @brief Removes every active condition reference in one family.
 *
 * Clearing a family does not change the stored objects themselves, their
 * amplitudes or the conditions in other families. The input reader must
 * independently clear or update any identifier-to-pointer references used
 * to resolve subsequent input-history operations.
 *
 * @param family Independent history family to clear.
 */
void ConditionManager::clear(ConditionFamily family) {
    conditions_[family].clear();
}

/**
 * @brief Returns the currently active condition pointers for one family.
 *
 * Model assembly reads this set directly rather than materializing duplicate
 * per-analysis condition lists. The reference remains tied to the manager's
 * stored family; insertions can invalidate iterators, and erasure invalidates
 * iterators to the erased element. The returned container must not be mutated
 * through this interface.
 *
 * @param family Family selected for assembly or membership inspection.
 * @return Constant reference to the live family set.
 */
const ConditionManager::Conditions& ConditionManager::get(ConditionFamily family) const {
    return conditions_[family];
}

} // namespace fem::bc
