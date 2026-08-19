/**
 * @file line2a.h
 * @brief Declares the two-node linear line element on [-1,1].
 */

#pragma once

#include "line.h"

namespace fem::model {

struct Line2A : public Line<2, MINUS_ONE_TO_ONE> {
    using Line<2, MINUS_ONE_TO_ONE>::Line;

    LinePtr copy() const override { return std::make_shared<Line2A>(*this); }

    StaticMatrix<2, 1> _shape_function(Precision r) const override;
    StaticMatrix<2, 1> shape_derivative(Precision r) const override;
    StaticMatrix<2, 1> shape_second_derivative(Precision r) const override;
};

} // namespace fem::model
