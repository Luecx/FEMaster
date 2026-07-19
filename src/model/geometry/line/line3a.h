/**
 * @file line3a.h
 * @brief Declares the three-node quadratic line element on [-1,1].
 *
 * `Line3A` uses quadratic Lagrange interpolation with two endpoint nodes and
 * one midpoint node on the symmetric natural interval.
 *
 * @see Line
 * @see LineInterface
 */

#pragma once

#include "line.h"

namespace fem::model {

/**
 * @brief Three-node quadratic line element on the natural interval `[-1,1]`.
 *
 * The first two nodes define the interval endpoints and the third node is the
 * midpoint. The element provides the quadratic shape functions and their first
 * and second natural-coordinate derivatives.
 */
struct Line3A : public Line<3, MINUS_ONE_TO_ONE> {
    using Line<3, MINUS_ONE_TO_ONE>::Line;

    // Quadratic Lagrange interpolation on the symmetric reference interval.
    // The first two nodes are the endpoints and the third node is the midpoint;
    // all derivatives follow this fixed local ordering.
    StaticMatrix<3,1> _shape_function(Precision r) const override;
    StaticMatrix<3,1> shape_derivative(Precision r) const override;
    StaticMatrix<3,1> shape_second_derivative(Precision r) const override;
};
}  // namespace fem::model
