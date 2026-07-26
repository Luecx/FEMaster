/**
 * @file rayleigh_damping.h
 * @brief Defines proportional Rayleigh damping for linear dynamic analyses.
 *
 * Rayleigh damping approximates the global viscous damping matrix as a linear
 * combination of the mass and stiffness matrices:
 *
 *     C = alpha M + beta K.
 *
 * The mass-proportional coefficient `alpha` primarily affects the lower
 * frequency range, while the stiffness-proportional coefficient `beta`
 * primarily affects the higher frequency range.
 *
 * The damping model stores only the two proportionality coefficients. The
 * global damping matrix is constructed from externally assembled mass and
 * stiffness matrices and is returned as a new sparse matrix.
 *
 * Matrix assembly, boundary-condition reduction and time integration remain
 * responsibilities of the surrounding load-case implementation.
 *
 * @see SparseMatrix
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem::loadcase::tools {

/**
 * @brief Proportional Rayleigh damping model.
 *
 * The model constructs the viscous damping matrix from the global mass and
 * stiffness matrices according to
 *
 *     C = alpha M + beta K.
 *
 * Both coefficients may be zero. In that case, the generated damping matrix
 * has the same dimensions as the supplied system matrices but contains no
 * entries.
 *
 * The class does not own or cache the mass, stiffness or damping matrices.
 * Every call to `build` constructs and returns a new sparse matrix.
 */
struct RayleighDamping {
    // Rayleigh damping coefficients
    Precision alpha = Precision(0);
    Precision beta  = Precision(0);

    // Construction
    RayleighDamping() = default;
    RayleighDamping(Precision alpha, Precision beta);

    // Construct the global damping matrix from the supplied mass and stiffness
    // matrices according to C = alpha M + beta K.
    [[nodiscard]] SparseMatrix build(const SparseMatrix& mass, const SparseMatrix& stiffness) const;
};

} // namespace fem::loadcase::tools