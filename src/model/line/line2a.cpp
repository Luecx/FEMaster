#include "line2a.h"

namespace fem::model {

Eigen::Matrix<Precision, 2, 1> Line2A::shape_function(Precision r) const {
    // Shape functions for a 2-node line element in the range [-1, 1]
    Eigen::Matrix<Precision, 2, 1> N;
    N(0) = 0.5 * (1 - r);
    N(1) = 0.5 * (1 + r);
    return N;
}

Eigen::Matrix<Precision, 2, 1> Line2A::shape_derivative(Precision /*r*/) const {
    // Shape function derivatives
    Eigen::Matrix<Precision, 2, 1> dN;
    dN(0) = -0.5;
    dN(1) = 0.5;
    return dN;
}

Eigen::Matrix<Precision, 2, 1> Line2A::shape_second_derivative(Precision /*r*/) const {
    // Second derivatives for a linear element are zero
    return Eigen::Matrix<Precision, 2, 1>::Zero();
}

}  // namespace fem::model
