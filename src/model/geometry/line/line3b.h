/**
 * @file line3b.h
 * @brief Declares the three-node quadratic line element on [0,1].
 *
 * `Line3B` uses quadratic Lagrange interpolation on `[0,1]` and applies the
 * corresponding first- and second-derivative scaling when mapping from the
 * symmetric reference interval.
 *
 * @see Line
 * @see LineInterface
 */

#pragma once

#include "line.h"

namespace fem::model {

/**
 * @brief Three-node quadratic line element on the natural interval `[0,1]`.
 *
 * The first two nodes define the interval endpoints and the third node is the
 * midpoint. The element provides quadratic interpolation and derivatives in
 * the non-symmetric natural coordinate system.
 */
struct Line3B : public Line<3, ZERO_TO_ONE> {
    using Line<3, ZERO_TO_ONE>::Line;

    // Quadratic Lagrange interpolation on [0,1]. The natural coordinate is
    // mapped to the symmetric interval for evaluation, and the first- and
    // second-derivative routines apply the corresponding chain-rule factors.
    StaticMatrix<3,1> _shape_function(Precision r) const override;
    StaticMatrix<3,1> shape_derivative(Precision r) const override;
    StaticMatrix<3,1> shape_second_derivative(Precision r) const override;
};
}  // namespace fem::model
