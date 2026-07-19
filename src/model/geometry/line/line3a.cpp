/**
 * @file line3a.cpp
 * @brief Implements quadratic Lagrange interpolation on the reference interval [-1,1].
 *
 * The element-specific routines provide the quadratic shape functions and
 * their first and second natural-coordinate derivatives for the common line
 * geometry implementation.
 *
 * @see Line3A
 */

#include "line3a.h"

namespace fem::model {

/**
 * Evaluates the three quadratic Lagrange shape functions on `[-1,1]`.
 *
 * The first two entries correspond to the interval endpoints and the third
 * entry corresponds to the midpoint node.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector of shape-function values.
 */
StaticMatrix<3,1> Line3A::_shape_function(Precision r) const {
    // Evaluate the endpoint and midpoint interpolation functions
    StaticMatrix<3,1> N;
    N(0) = 0.5 * r * (r - 1);  // N1
    N(1) = 0.5 * r * (r + 1);  // N2
    N(2) = (1 - r) * (1 + r);  // N3 (midpoint node)
    return N;
}

/**
 * Evaluates the first natural-coordinate derivatives of the quadratic shape
 * functions.
 *
 * @param r Natural line coordinate.
 * @return Fixed-size vector containing `dN/dr` for all three nodes.
 */
StaticMatrix<3,1> Line3A::shape_derivative(Precision r) const {
    // Evaluate the derivatives of the endpoint and midpoint functions
    StaticMatrix<3,1> dN;
    dN(0) = r - 0.5;   // dN1/dr
    dN(1) = r + 0.5;   // dN2/dr
    dN(2) = -2.0 * r;  // dN3/dr
    return dN;
}

/**
 * Evaluates the constant second derivatives of the quadratic shape functions.
 *
 * @param r Natural line coordinate, unused because the second derivatives are constant.
 * @return Fixed-size vector containing `d²N/dr²` for all three nodes.
 */
StaticMatrix<3,1> Line3A::shape_second_derivative(Precision /*r*/) const {
    // Store the constant second derivatives of the quadratic functions
    StaticMatrix<3,1> ddN;
    ddN(0) = 1.0;   // d^2N1/dr^2
    ddN(1) = 1.0;   // d^2N2/dr^2
    ddN(2) = -2.0;  // d^2N3/dr^2
    return ddN;
}
}  // namespace fem::model
