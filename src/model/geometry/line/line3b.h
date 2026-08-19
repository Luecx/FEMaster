/**
 * @file line3b.h
 * @brief Declares the three-node quadratic line element on [0,1].
 */

#pragma once

#include "line.h"

namespace fem::model {

struct Line3B : public Line<3, ZERO_TO_ONE> {
    using Line<3, ZERO_TO_ONE>::Line;

    LinePtr copy() const override { return std::make_shared<Line3B>(*this); }

    StaticMatrix<3, 1> _shape_function(Precision r) const override;
    StaticMatrix<3, 1> shape_derivative(Precision r) const override;
    StaticMatrix<3, 1> shape_second_derivative(Precision r) const override;
};

} // namespace fem::model
