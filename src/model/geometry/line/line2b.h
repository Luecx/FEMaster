/**
 * @file line2b.h
 * @brief Declares the two-node linear line element on the reference interval [0,1].
 *
 * `Line2B` uses linear Lagrange interpolation on `[0,1]`. Its implementation
 * maps this interval to the symmetric reference interval used by the shape
 * functions while preserving the correct derivative scaling.
 *
 * @see Line
 * @see LineInterface
 */

#pragma once

#include "line.h"

namespace fem::model {

/**
 * @brief Two-node linear line element on the natural interval `[0,1]`.
 *
 * The two nodes occupy the endpoints of the reference interval. The shape
 * functions are linear, their first derivatives include the interval mapping
 * factor, and all second derivatives vanish.
 */
struct Line2B : public Line<2, ZERO_TO_ONE> {
    using Line<2, ZERO_TO_ONE>::Line;

    // Linear Lagrange interpolation on [0,1]. The implementation maps the
    // natural coordinate to [-1,1] and applies the corresponding chain-rule
    // factor to the first derivatives before they are used by line geometry.
    StaticMatrix<2,1> _shape_function(Precision r) const override;
    StaticMatrix<2,1> shape_derivative(Precision r) const override;
    StaticMatrix<2,1> shape_second_derivative(Precision r) const override;
};
}  // namespace fem::model
