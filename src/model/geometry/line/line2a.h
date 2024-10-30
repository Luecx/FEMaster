#pragma once

#include "line_interface.h"

namespace fem::model {

/**
 * @brief Class representing a 2-node line element with coordinate range [-1, 1].
 */
struct Line2A : public LineInterface<2, MINUS_ONE_TO_ONE> {
    using LineInterface<2, MINUS_ONE_TO_ONE>::LineInterface;

    StaticMatrix<2,1> shape_function(Precision r) const override;
    StaticMatrix<2,1> shape_derivative(Precision r) const override;
    StaticMatrix<2,1> shape_second_derivative(Precision r) const override;
};

}  // namespace fem::model
