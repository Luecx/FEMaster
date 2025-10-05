/**
 * @file strain.h
 * @brief Declares the strain tensor utility in Voigt form.
 *
 * Provides transformation helpers to rotate strain vectors between coordinate
 * systems.
 *
 * @see src/material/strain.cpp
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
 * @struct Strain
 * @brief Represents strain in Voigt notation and offers rotation helpers.
 */
struct Strain : Vec6 {
    /**
     * @brief Rotates the strain vector into a different basis.
     *
     * @param basis Target coordinate basis.
     * @return Strain Rotated strain vector.
     */
    [[nodiscard]] Strain transform(const cos::Basis& basis) const;

    /**
     * @brief Computes the 6x6 transformation matrix for a given basis.
     *
     * @param basis Target coordinate basis.
     * @return Mat6 Transformation matrix in Voigt notation.
     */
    static Mat6 get_transformation_matrix(const cos::Basis& basis);
};

using Strains = std::vector<Strain>; ///< Collection alias for strain vectors.

} // namespace fem
