/**
 * @file feature.h
 * @brief Base interface for non-element features contributing to global operators.
 *
 * Features add model behavior that is not represented by regular finite
 * elements, for example concentrated point masses or springs.
 *
 * @see src/feature/point_mass.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "../core/types_eig.h"

namespace fem {
namespace feature {
/**
 * @struct Feature
 * @brief Abstract base interface for non-element matrix contributions.
 *
 * Implementations assemble directly into the global stiffness or mass triplet
 * list, but only for already active system DOFs.
 */
struct Feature {
    using Ptr = std::shared_ptr<Feature>; ///< Shared pointer alias for feature ownership.

    /**
     * @brief Defaulted virtual destructor for polymorphic deletion.
     */
    virtual ~Feature() = default;

    /**
     * @brief Adds this feature's stiffness terms to a global triplet list.
     *
     * @param indices Active global DOF ids. Negative entries are inactive.
     * @param out Triplet list receiving stiffness entries.
     */
    virtual void assemble_stiffness(const SystemDofIds& indices, TripletList& out) const = 0;

    /**
     * @brief Adds this feature's mass terms to a global triplet list.
     *
     * @param indices Active global DOF ids. Negative entries are inactive.
     * @param out Triplet list receiving mass entries.
     */
    virtual void assemble_mass(const SystemDofIds& indices, TripletList& out) const = 0;
};
} // namespace feature
} // namespace fem
