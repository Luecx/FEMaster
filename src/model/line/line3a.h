#pragma once

#include "line_interface.h"

namespace fem::model {

/**
 * @brief Class representing a 3-node quadratic line element with coordinate range [-1, 1].
 */
struct Line3A : public LineInterface<3, MINUS_ONE_TO_ONE> {
    using LineInterface<3, MINUS_ONE_TO_ONE>::LineInterface;

    Eigen::Matrix<Precision, 3, 1> shape_function(Precision r) const override;
    Eigen::Matrix<Precision, 3, 1> shape_derivative(Precision r) const override;
    Eigen::Matrix<Precision, 3, 1> shape_second_derivative(Precision r) const override;
};

}  // namespace fem::model
