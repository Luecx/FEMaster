#include "line3a.h"

namespace fem::model {

StaticMatrix<3,1> Line3A::shape_function(Precision r) const {
    // Shape functions for a 3-node quadratic line element in the range [-1, 1]
    StaticMatrix<3,1> N;
    N(0) = 0.5 * r * (r - 1);  // N1
    N(1) = 0.5 * r * (r + 1);  // N2
    N(2) = (1 - r) * (1 + r);  // N3 (midpoint node)
    return N;
}

StaticMatrix<3,1> Line3A::shape_derivative(Precision r) const {
    // Shape function derivatives for a 3-node quadratic line element
    StaticMatrix<3,1> dN;
    dN(0) = r - 0.5;   // dN1/dr
    dN(1) = r + 0.5;   // dN2/dr
    dN(2) = -2.0 * r;  // dN3/dr
    return dN;
}

StaticMatrix<3,1> Line3A::shape_second_derivative(Precision /*r*/) const {
    // Second derivatives for a 3-node quadratic element
    StaticMatrix<3,1> ddN;
    ddN(0) = 1.0;   // d^2N1/dr^2
    ddN(1) = 1.0;   // d^2N2/dr^2
    ddN(2) = -2.0;  // d^2N3/dr^2
    return ddN;
}

}  // namespace fem::model
