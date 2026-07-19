/**
 * @file line2a.h
 * @brief Declares the two-node linear line element on the reference interval [-1,1].
 *
 * `Line2A` uses linear Lagrange interpolation with the natural coordinate
 * range `[-1,1]` and provides the element-specific shape functions required by
 * the common line geometry implementation.
 *
 * @see Line
 * @see LineInterface
 */

#pragma once

#include "line.h"

namespace fem::model {

/**
 * @brief Two-node linear line element on the natural interval `[-1,1]`.
 *
 * The two nodes occupy the endpoints of the reference interval. The shape
 * functions are linear, their first derivatives are constant, and all second
 * derivatives vanish.
 */
struct Line2A : public Line<2, MINUS_ONE_TO_ONE> {
    using Line<2, MINUS_ONE_TO_ONE>::Line;

    // Linear Lagrange interpolation on the symmetric reference interval. The
    // underscore overload supplies fixed-size values, while the derivative
    // routines provide the constant natural-coordinate tangent coefficients
    // required by the common line geometry implementation.
    StaticMatrix<2,1> _shape_function(Precision r) const override;
    StaticMatrix<2,1> shape_derivative(Precision r) const override;
    StaticMatrix<2,1> shape_second_derivative(Precision r) const override;
};
}  // namespace fem::model
