#include "line3b.h"

namespace fem::model {

Eigen::Matrix<Precision, 3, 1> Line3B::shape_function(Precision r) const {
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    Eigen::Matrix<Precision, 3, 1> N;
    N(0) = 0.5 * r * (r - 1);  // N1
    N(1) = 0.5 * r * (r + 1);  // N2
    N(2) = (1 - r) * (1 + r);  // N3 (midpoint node)
    return N;
}

Eigen::Matrix<Precision, 3, 1> Line3B::shape_derivative(Precision r) const {
    r = 2.0 * r - 1.0;  // Map from [0, 1] to [-1, 1]
    Eigen::Matrix<Precision, 3, 1> dN;
    dN(0) = r - 0.5;   // dN1/dr
    dN(1) = r + 0.5;   // dN2/dr
    dN(2) = -2.0 * r;  // dN3/dr

    dN *= 2.0;  // Scale derivatives by the mapping factor

    return dN;
}

Eigen::Matrix<Precision, 3, 1> Line3B::shape_second_derivative(Precision /*r*/) const {
    Eigen::Matrix<Precision, 3, 1> ddN;
    ddN(0) = 1.0;   // d^2N1/dr^2
    ddN(1) = 1.0;   // d^2N2/dr^2
    ddN(2) = -2.0;  // d^2N3/dr^2

    ddN *= 4.0;  // Scale second derivatives by the mapping factor

    return ddN;
}

}  // namespace fem::model
