#pragma once

#include "line_interface.h"

namespace fem::model {

/**
 * @brief Class representing a 3-node quadratic line element with coordinate range [0, 1].
 */
struct Line3B : public LineInterface<3, ZERO_TO_ONE> {
    using LineInterface<3, ZERO_TO_ONE>::LineInterface;

    StaticMatrix<3,1> shape_function(Precision r) const override;
    StaticMatrix<3,1> shape_derivative(Precision r) const override;
    StaticMatrix<3,1> shape_second_derivative(Precision r) const override;
};

}  // namespace fem::model
