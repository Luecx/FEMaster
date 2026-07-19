/**
 * @file line2b.cpp
 * @brief Implements linear Lagrange interpolation on the reference interval [0,1].
 *
 * The shape functions are evaluated through the affine mapping from `[0,1]`
 * to the symmetric interval `[-1,1]`; derivative routines include the required
 * chain-rule scaling.
 *
 * @see Line2B
 */

#include "line2b.h"

namespace fem::model {

/**
 * Evaluates the two linear shape functions on `[0,1]`.
 *
 * The natural coordinate is first mapped to `[-1,1]`, where the standard
 * endpoint Lagrange functions are evaluated.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector of shape-function values.
 */
StaticMatrix<2,1> Line2B::_shape_function(Precision r) const {
    // Map the natural coordinate to the symmetric reference interval
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    StaticMatrix<2,1> N;
    N(0) = 0.5 * (1 - r);
    N(1) = 0.5 * (1 + r);
    return N;
}

/**
 * Evaluates the first derivatives of the linear shape functions on `[0,1]`.
 *
 * The derivatives are scaled by two because the symmetric coordinate is
 * `r_sym = 2 r - 1`.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector containing `dN/dr` for both nodes.
 */
StaticMatrix<2,1> Line2B::shape_derivative(Precision r) const {
    // Map the coordinate before evaluating the symmetric-interval derivatives
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    StaticMatrix<2,1> dN;
    dN(0) = -0.5;
    dN(1) = 0.5;

    // Apply the first-derivative chain-rule factor of the affine mapping
    dN *= 2.0;

    return dN;
}

/**
 * Evaluates the second derivatives of the linear shape functions.
 *
 * Linear interpolation has no quadratic contribution, so the second
 * derivatives remain zero after the interval mapping.
 *
 * @param r Natural line coordinate, unused.
 * @return Zero vector containing the second derivatives.
 */
StaticMatrix<2,1> Line2B::shape_second_derivative(Precision /*r*/) const {
    // Linear interpolation has identically zero second derivatives
    return StaticMatrix<2,1>::Zero();
}
}  // namespace fem::model
