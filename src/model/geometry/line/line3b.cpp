/**
 * @file line3b.cpp
 * @brief Implements quadratic Lagrange interpolation on the reference interval [0,1].
 *
 * The shape functions are evaluated through the affine mapping from `[0,1]`
 * to `[-1,1]`. The first and second derivative routines include the
 * corresponding chain-rule scaling factors.
 *
 * @see Line3B
 */

#include "line3b.h"

namespace fem::model {

/**
 * Evaluates the three quadratic Lagrange shape functions on `[0,1]`.
 *
 * The natural coordinate is mapped to `[-1,1]` before evaluating the standard
 * endpoint and midpoint interpolation functions.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector of shape-function values.
 */
StaticMatrix<3,1> Line3B::_shape_function(Precision r) const {
    // Map the natural coordinate to the symmetric reference interval
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    StaticMatrix<3,1> N;
    N(0) = 0.5 * r * (r - 1);  // N1
    N(1) = 0.5 * r * (r + 1);  // N2
    N(2) = (1 - r) * (1 + r);  // N3 (midpoint node)
    return N;
}

/**
 * Evaluates the first natural-coordinate derivatives of the quadratic shape
 * functions on `[0,1]`.
 *
 * The derivatives are multiplied by two for the affine coordinate mapping.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector containing `dN/dr` for all three nodes.
 */
StaticMatrix<3,1> Line3B::shape_derivative(Precision r) const {
    // Map the coordinate before evaluating the symmetric-interval derivatives
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    StaticMatrix<3,1> dN;
    dN(0) = r - 0.5;   // dN1/dr
    dN(1) = r + 0.5;   // dN2/dr
    dN(2) = -2.0 * r;  // dN3/dr

    // Apply the first-derivative chain-rule factor of the affine mapping
    dN *= 2.0;

    return dN;
}

/**
 * Evaluates the second natural-coordinate derivatives of the quadratic shape
 * functions on `[0,1]`.
 *
 * The second derivatives are multiplied by four because two derivatives are
 * taken through the affine coordinate mapping.
 *
 * @param r Natural line coordinate, unused because the second derivatives are constant.
 * @return Fixed-size vector containing `d²N/dr²` for all three nodes.
 */
StaticMatrix<3,1> Line3B::shape_second_derivative(Precision /*r*/) const {
    // Store the constant second derivatives on the symmetric interval
    StaticMatrix<3,1> ddN;
    ddN(0) = 1.0;   // d^2N1/dr^2
    ddN(1) = 1.0;   // d^2N2/dr^2
    ddN(2) = -2.0;  // d^2N3/dr^2

    // Apply the second-derivative chain-rule factor of the affine mapping
    ddN *= 4.0;

    return ddN;
}
}  // namespace fem::model
