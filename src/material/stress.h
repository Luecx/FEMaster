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
 *
 * Voigt ordering:
 *
 *     [s11, s22, s33, s23, s13, s12]^T
 *
 * The transformation is interpreted as a passive basis change:
 *
 *     stress_to = T(from_basis -> to_basis) * stress_from
 */
struct Stress : Vec6 {
    /**
     * @brief Expresses this stress vector in another basis.
     *
     * @param from_basis Basis in which this stress vector is currently expressed.
     * @param to_basis Target basis in which the returned stress vector is expressed.
     * @return Stress Stress vector expressed in @p to_basis.
     */
    [[nodiscard]] Stress transform(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    ) const;

    /**
     * @brief Computes the 6x6 stress transformation matrix between two bases.
     *
     * @param from_basis Source basis.
     * @param to_basis Target basis.
     * @return Mat6 Transformation matrix in Voigt notation.
     */
    static Mat6 get_transformation_matrix(
        const cos::Basis& from_basis,
        const cos::Basis& to_basis
    );
};

using Stresses = std::vector<Stress>; ///< Collection alias for stress vectors.

} // namespace fem