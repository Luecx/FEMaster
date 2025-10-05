/**
 * @file stress.h
 * @brief Declares the stress tensor utility in Voigt form.
 *
 * Provides transformation helpers to rotate stress vectors between coordinate
 * systems.
 *
 * @see src/material/stress.cpp
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"

#include <vector>

namespace fem {

/**
 * @struct Stress
 * @brief Represents stress in Voigt notation and offers rotation helpers.
 */
struct Stress : Vec6 {
    /**
     * @brief Rotates the stress vector into a different basis.
     *
     * @param basis Target coordinate basis.
     * @return Stress Rotated stress vector.
     */
    [[nodiscard]] Stress transform(const cos::Basis& basis) const;

    /**
     * @brief Computes the 6x6 transformation matrix for a given basis.
     *
     * @param basis Target coordinate basis.
     * @return Mat6 Transformation matrix in Voigt notation.
     */
    static Mat6 get_transformation_matrix(const cos::Basis& basis);
};

using Stresses = std::vector<Stress>; ///< Collection alias for stress vectors.

} // namespace fem
