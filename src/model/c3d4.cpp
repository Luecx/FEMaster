//
// Created by Luecx on 12.06.2023.
//

#include "c3d4.h"

namespace fem {
namespace model {

C3D4::C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

StaticMatrix<4, 1> C3D4::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<4, 1>  res {};

    res(0) = r;
    res(1) = s;
    res(2) = 1 - r - s - t;
    res(3) = t;

    return res;
}

StaticMatrix<4, 3> C3D4::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<4, 3> local_shape_derivative {};
    local_shape_derivative.setZero();

    local_shape_derivative(0,0) = 1;

    local_shape_derivative(1,1) = 1;

    local_shape_derivative(2,0) = -1;
    local_shape_derivative(2,1) = -1;
    local_shape_derivative(2,2) = -1;

    local_shape_derivative(3,2) = 1;

    return local_shape_derivative;
}

StaticMatrix<4, 3> C3D4::node_coords_local() {
    StaticMatrix<4, 3> res {};
    res.setZero();
    res(0, 0) = 1;
    res(0, 1) = 0;
    res(0, 2) = 0;

    res(1, 0) = 0;
    res(1, 1) = 1;
    res(1, 2) = 0;

    res(2, 0) = 0;
    res(2, 1) = 0;
    res(2, 2) = 0;

    res(3, 0) = 0;
    res(3, 1) = 0;
    res(3, 2) = 1;
    return res;
}

}    // namespace model
}    // namespace fem