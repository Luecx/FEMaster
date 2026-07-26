/**
 * @file rayleigh_damping.cpp
 * @brief Implements proportional Rayleigh damping.
 *
 * The implementation constructs a sparse viscous damping matrix from the
 * supplied mass and stiffness matrices according to
 *
 *     C = alpha M + beta K.
 *
 * Dedicated branches handle vanishing coefficients to avoid unnecessary
 * sparse-matrix operations. When both coefficients are zero, an empty sparse
 * matrix with the dimensions of the supplied system matrices is returned.
 *
 * @see RayleighDamping
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#include "rayleigh_damping.h"

namespace fem::loadcase::tools {

RayleighDamping::RayleighDamping(Precision alpha, Precision beta)
    : alpha(alpha),
      beta (beta) {}

/**
 * Constructs the proportional viscous damping matrix.
 *
 * The supplied matrices are combined according to
 *
 *     C = alpha M + beta K.
 *
 * Separate branches are used for the cases in which one or both damping
 * coefficients vanish:
 *
 * - if both coefficients are zero, an empty matrix is returned,
 * - if `alpha` is zero, only the stiffness-proportional contribution is built,
 * - if `beta` is zero, only the mass-proportional contribution is built,
 * - otherwise, both sparse contributions are assembled.
 *
 * The mass and stiffness matrices are expected to describe the same equation
 * system and therefore to have matching dimensions. Their structural
 * compatibility is assumed to have been established by the calling load case.
 *
 * @param mass Global or reduced mass matrix.
 * @param stiffness Global or reduced stiffness matrix.
 *
 * @return Newly constructed sparse damping matrix.
 */
SparseMatrix RayleighDamping::build(const SparseMatrix& mass, const SparseMatrix& stiffness) const {
    // Return an empty matrix when damping is completely disabled
    if (alpha == Precision(0) && beta == Precision(0)) {
        return SparseMatrix(mass.rows(), mass.cols());
    }

    // Avoid evaluating the mass contribution when only stiffness-proportional
    // damping is active
    if (alpha == Precision(0)) {
        return beta * stiffness;
    }

    // Avoid evaluating the stiffness contribution when only mass-proportional
    // damping is active
    if (beta == Precision(0)) {
        return alpha * mass;
    }

    // Assemble both proportional damping contributions
    SparseMatrix damping = alpha * mass;
    damping             += beta * stiffness;

    damping.makeCompressed();

    return damping;
}

} // namespace fem::loadcase::tools