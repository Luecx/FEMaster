/**
 * @file line2a.cpp
 * @brief Implements linear Lagrange interpolation on the reference interval [-1,1].
 *
 * The element-specific routines provide the shape functions and their natural
 * coordinate derivatives used by the generic line geometry algorithms.
 *
 * @see Line2A
 */

#include "line2a.h"

namespace fem::model {

/**
 * Evaluates the two linear shape functions on `[-1,1]`.
 *
 * The first and second entries correspond to the nodes at `r=-1` and `r=1`.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector of shape-function values.
 */
StaticMatrix<2,1> Line2A::_shape_function(Precision r) const {
    // Evaluate the endpoint interpolation functions on the symmetric interval
    StaticMatrix<2,1> N;
    N(0) = 0.5 * (1 - r);
    N(1) = 0.5 * (1 + r);
    return N;
}

/**
 * Evaluates the constant first derivatives of the linear shape functions.
 *
 * @param r Natural line coordinate, unused because the derivatives are constant.
 * @return Fixed-size vector containing `dN/dr` for both nodes.
 */
StaticMatrix<2,1> Line2A::shape_derivative(Precision /*r*/) const {
    // Store the constant natural-coordinate derivatives
    StaticMatrix<2,1> dN;
    dN(0) = -0.5;
    dN(1) = 0.5;
    return dN;
}

/**
 * Evaluates the second derivatives of the linear shape functions.
 *
 * Linear interpolation has no quadratic contribution, so all second
 * derivatives vanish.
 *
 * @param r Natural line coordinate, unused.
 * @return Zero vector containing the second derivatives.
 */
StaticMatrix<2,1> Line2A::shape_second_derivative(Precision /*r*/) const {
    // Linear interpolation has identically zero second derivatives
    return StaticMatrix<2,1>::Zero();
}
}  // namespace fem::model
