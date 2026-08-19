/**
 * @file line3a.h
 * @brief Declares the three-node quadratic line element on [-1,1].
 */

#pragma once

#include "line.h"

namespace fem::model {

struct Line3A : public Line<3, MINUS_ONE_TO_ONE> {
    using Line<3, MINUS_ONE_TO_ONE>::Line;

    LinePtr copy() const override { return std::make_shared<Line3A>(*this); }

    StaticMatrix<3, 1> _shape_function(Precision r) const override;
    StaticMatrix<3, 1> shape_derivative(Precision r) const override;
    StaticMatrix<3, 1> shape_second_derivative(Precision r) const override;
};

} // namespace fem::model
