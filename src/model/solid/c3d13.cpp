// Created by Finn on 19.10.2024

#include "c3d13.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface3.h"

namespace fem {
namespace model {

C3D13::C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids)
    : SolidElement(p_elem_id, p_node_ids) {}

const quadrature::Quadrature& C3D13::integration_scheme() {
    const static quadrature::Quadrature quad {quadrature::DOMAIN_ISO_PYRAMID, quadrature::ORDER_QUARTIC};
    return quad;
}
SurfacePtr C3D13::surface(ID surface_id) {
    // TODO: implement this
    (void) surface_id;
    return nullptr;
}

StaticMatrix<13, 1> C3D13::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<13, 1> res {};

    // Define shape functions as provided
    res(0, 0) = 0.25 * (r + s - 1) * ((1 + r) * (1 + s) - t + (r * s * t) / (1 - t));
    res(1, 0) = 0.25 * (-r + s - 1) * ((1 - r) * (1 + s) - t - (r * s * t) / (1 - t));
    res(2, 0) = 0.25 * (-r - s - 1) * ((1 - r) * (1 - s) - t + (r * s * t) / (1 - t));
    res(3, 0) = 0.25 * (r - s - 1) * ((1 + r) * (1 - s) - t - (r * s * t) / (1 - t));
    res(4, 0) = t * (2 * t - 1);
    res(5, 0) = ((1 + r - t) * (1 - r - t) * (1 + s - t)) / (2 * (1 - t));
    res(6, 0) = ((1 + s - t) * (1 - s - t) * (1 - r - t)) / (2 * (1 - t));
    res(7, 0) = ((1 + r - t) * (1 - r - t) * (1 - s - t)) / (2 * (1 - t));
    res(8, 0) = ((1 + s - t) * (1 - s - t) * (1 + r - t)) / (2 * (1 - t));
    res(9, 0) = t * (1 + r - t) * (1 + s - t) / (1 - t);
    res(10, 0) = t * (1 - r - t) * (1 + s - t) / (1 - t);
    res(11, 0) = t * (1 - r - t) * (1 - s - t) / (1 - t);
    res(12, 0) = t * (1 + r - t) * (1 - s - t) / (1 - t);

    return res;
}

StaticMatrix<13, 3> C3D13::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<13, 3> der {};

    // Define shape function derivatives as provided
    der(0, 0) = 0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (r + 1) * (s + 1)
                + (0.25 * r + 0.25 * s - 0.25) * (s * t / (1 - t) + s + 1);
    der(0, 1) = 0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (r + 1) * (s + 1)
                + (0.25 * r + 0.25 * s - 0.25) * (r * t / (1 - t) + r + 1);
    der(0, 2) = (0.25 * r + 0.25 * s - 0.25) * (r * s * t / (1 - t) / (1 - t) + r * s / (1 - t) - 1);

    der(1, 0) = 0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (s + 1)
                + (-0.25 * r + 0.25 * s - 0.25) * (-s * t / (1 - t) - s - 1);
    der(1, 1) = -0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (1 - r) * (s + 1)
                + (-0.25 * r + 0.25 * s - 0.25) * (-r * t / (1 - t) - r + 1);
    der(1, 2) = (-0.25 * r + 0.25 * s - 0.25) * (-r * s * t / (1 - t) / (1 - t) - r * s / (1 - t) - 1);

    der(2, 0) = -0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (1 - s)
                + (-0.25 * r - 0.25 * s - 0.25) * (s * t / (1 - t) + s - 1);
    der(2, 1) = -0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - r) * (1 - s)
                + (-0.25 * r - 0.25 * s - 0.25) * (r * t / (1 - t) + r - 1);
    der(2, 2) = (-0.25 * r - 0.25 * s - 0.25) * (r * s * t / (1 - t) / (1 - t) + r * s / (1 - t) - 1);

    der(3, 0) = -0.25 * r * s * t / (1 - t) - 0.25 * t + 0.25 * (1 - s) * (r + 1)
                + (0.25 * r - 0.25 * s - 0.25) * (-s * t / (1 - t) - s + 1);
    der(3, 1) = 0.25 * r * s * t / (1 - t) + 0.25 * t - 0.25 * (1 - s) * (r + 1)
                + (0.25 * r - 0.25 * s - 0.25) * (-r * t / (1 - t) - r - 1);
    der(3, 2) = (0.25 * r - 0.25 * s - 0.25) * (-r * s * t / (1 - t) / (1 - t) - r * s / (1 - t) - 1);

    der(4, 0) = 0;
    der(4, 1) = 0;
    der(4, 2) = 4 * t - 1;

    der(5, 0) = (-r - t + 1) * (s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(5, 1) = (-r - t + 1) * (r - t + 1) / (2 - 2 * t);
    der(5, 2) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (r - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (r - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(6, 0) = -(-s - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(6, 1) = (-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(6, 2) = -(-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (-r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (-s - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (-s - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(7, 0) = (-r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (-s - t + 1) / (2 - 2 * t);
    der(7, 1) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t);
    der(7, 2) = -(-r - t + 1) * (r - t + 1) / (2 - 2 * t) - (-r - t + 1) * (-s - t + 1) / (2 - 2 * t)
                - (r - t + 1) * (-s - t + 1) / (2 - 2 * t)
                + 2 * (-r - t + 1) * (r - t + 1) * (-s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(8, 0) = (-s - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(8, 1) = (r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t);
    der(8, 2) = -(r - t + 1) * (-s - t + 1) / (2 - 2 * t) - (r - t + 1) * (s - t + 1) / (2 - 2 * t)
                - (-s - t + 1) * (s - t + 1) / (2 - 2 * t)
                + 2 * (r - t + 1) * (-s - t + 1) * (s - t + 1) / ((2 - 2 * t) * (2 - 2 * t));

    der(9, 0) = t * (s - t + 1) / (1 - t);
    der(9, 1) = t * (r - t + 1) / (1 - t);
    der(9, 2) = -t * (r - t + 1) / (1 - t) - t * (s - t + 1) / (1 - t)
                + t * (r - t + 1) * (s - t + 1) / ((1 - t) * (1 - t)) + (r - t + 1) * (s - t + 1) / (1 - t);

    der(10, 0) = -t * (s - t + 1) / (1 - t);
    der(10, 1) = t * (-r - t + 1) / (1 - t);
    der(10, 2) = -t * (-r - t + 1) / (1 - t) - t * (s - t + 1) / (1 - t)
                 + t * (-r - t + 1) * (s - t + 1) / ((1 - t) * (1 - t)) + (-r - t + 1) * (s - t + 1) / (1 - t);

    der(11, 0) = -t * (-s - t + 1) / (1 - t);
    der(11, 1) = -t * (-r - t + 1) / (1 - t);
    der(11, 2) = -t * (-r - t + 1) / (1 - t) - t * (-s - t + 1) / (1 - t)
                 + t * (-r - t + 1) * (-s - t + 1) / ((1 - t) * (1 - t)) + (-r - t + 1) * (-s - t + 1) / (1 - t);

    der(12, 0) = t * (-s - t + 1) / (1 - t);
    der(12, 1) = -t * (r - t + 1) / (1 - t);
    der(12, 2) = -t * (r - t + 1) / (1 - t) - t * (-s - t + 1) / (1 - t)
                 + t * (r - t + 1) * (-s - t + 1) / ((1 - t) * (1 - t)) + (r - t + 1) * (-s - t + 1) / (1 - t);

    return der;
}

StaticMatrix<13, 3> C3D13::node_coords_local() {
    StaticMatrix<13, 3> res {};
    res <<   1  ,  1  , 0,
            -1  ,  1  , 0,
            -1  , -1  , 0,
             1  , -1  , 0,
             0  ,  0  , 1,
             0  ,  1  , 0,
            -1  ,  0  , 0,
             0  , -1  , 0,
             1  ,  0  , 0,
             0.5,  0.5, 0.5,
            -0.5,  0.5, 0.5,
            -0.5, -0.5, 0.5,
             0.5, -0.5, 0.5;

    return res;
}

}    // namespace model
}    // namespace fem
