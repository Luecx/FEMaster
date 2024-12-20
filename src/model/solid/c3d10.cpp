//
// Created by Luecx on 12.06.2023.
//

#include "c3d10.h"
#include "../geometry/surface/surface6.h"

namespace fem {
namespace model {

C3D10::C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

StaticMatrix<10, 1> C3D10::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<10, 1> res {};

    // Vertex nodes
    res(0) = (1 - r - s - t) * (2*(1 - r - s - t) - 1);
    res(1) = r * (2*r - 1);
    res(2) = s * (2*s - 1);
    res(3) = t * (2*t - 1);

    // Mid-edge nodes
    res(4) = 4 * r * (1 - r - s - t);
    res(5) = 4 * r * s;
    res(6) = 4 * s * (1 - r - s - t);
    res(7) = 4 * t * (1 - r - s - t);
    res(8) = 4 * r * t;
    res(9) = 4 * s * t;

    return res;
}


StaticMatrix<10, 3> C3D10::shape_derivative(Precision r, Precision s, Precision t) {
    StaticMatrix<10, 3> res {};
    res.setZero();

    Precision u = 1 - r - s - t;

    //    res(0) = r * (2*r - 1);
    //    res(1) = s * (2*s - 1);
    //    res(2) = (1 - r - s - t) * (2*(1 - r - s - t) - 1);
    //    res(3) = t * (2*t - 1);
    //
    //    // Mid-edge nodes
    //    res(4) = 4 * r * s;                    // Node 5 (between 1-2)
    //    res(5) = 4 * s * t;                    // Node 6 (flipped, now between 2-3)
    //    res(6) = 4 * r * t;                    // Node 7 (flipped, now between 1-3)
    //    res(7) = 4 * r * (1 - r - s - t);      // Node 8 (between 1-4)
    //    res(8) = 4 * s * (1 - r - s - t);      // Node 9 (between 2-4)
    //    res(9) = 4 * t * (1 - r - s - t);      // Node 10 (between 3-4)

//    res(0) = (1 - r - s - t) * (2*(1 - r - s - t) - 1);
//    res(1) = r * (2*r - 1);
//    res(2) = s * (2*s - 1);
//    res(3) = t * (2*t - 1);
//
//    // Mid-edge nodes
//    res(4) = 4 * r * (1 - r - s - t);
//    res(5) = 4 * r * s;
//    res(6) = 4 * s * (1 - r - s - t);
//    res(7) = 4 * t * (1 - r - s - t);
//    res(8) = 4 * r * t;
//    res(9) = 4 * s * t;

    // Derivatives with respect to r
    res(0, 0) = -4*u + 1;
    res(1, 0) = 4*r - 1;
    res(2, 0) = 0;
    res(3, 0) = 0;
    res(4, 0) = 4*u - 4*r;
    res(5, 0) = 4*s;
    res(6, 0) = -4*s;
    res(7, 0) = -4*t;
    res(8, 0) = 4*t;
    res(9, 0) = 0;

    // Derivatives with respect to s
    res(0, 1) = -4*u + 1;
    res(1, 1) = 0;
    res(2, 1) = 4*s - 1;
    res(3, 1) = 0;
    res(4, 1) = -4*r;
    res(5, 1) = 4*r;
    res(6, 1) = 4*u - 4*s;
    res(7, 1) = -4*t;
    res(8, 1) = 0;
    res(9, 1) = 4*t;

    // Derivatives with respect to t
    res(0, 2) = -4*u + 1;
    res(1, 2) = 0;
    res(2, 2) = 0;
    res(3, 2) = 4*t - 1;
    res(4, 2) = -4*r;
    res(5, 2) = 0;
    res(6, 2) = -4*s;
    res(7, 2) = 4*u - 4*t;
    res(8, 2) = 4*r;
    res(9, 2) = 4*s;

    return res;
}


StaticMatrix<10, 3> C3D10::node_coords_local() {
    StaticMatrix<10, 3> res {};
    res.setZero();

    // Vertex nodes
    res(0, 0) = 0;
    res(0, 1) = 0;
    res(0, 2) = 0;
    res(1, 0) = 1;
    res(1, 1) = 0;
    res(1, 2) = 0;
    res(2, 0) = 0;
    res(2, 1) = 1;
    res(2, 2) = 0;
    res(3, 0) = 0;
    res(3, 1) = 0;
    res(3, 2) = 1;

    // Mid-edge nodes
    res(4, 0) = 0.5;
    res(4, 1) = 0.0;
    res(4, 2) = 0;
    res(5, 0) = 0.5;
    res(5, 1) = 0.5;
    res(5, 2) = 0.0;
    res(6, 0) = 0.0;
    res(6, 1) = 0.5;
    res(6, 2) = 0.0;
    res(7, 0) = 0.0;
    res(7, 1) = 0.0;
    res(7, 2) = 0.5;
    res(8, 0) = 0.5;
    res(8, 1) = 0.0;
    res(8, 2) = 0.5;
    res(9, 0) = 0.0;
    res(9, 1) = 0.5;
    res(9, 2) = 0.5;

    return res;
}
SurfacePtr C3D10::surface(ID surface_id) {
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[0], node_ids[1], node_ids[2], node_ids[4], node_ids[5], node_ids[6]});
        case 2:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[0], node_ids[3], node_ids[1], node_ids[7], node_ids[8], node_ids[4]});
        case 3:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[1], node_ids[3], node_ids[2], node_ids[8], node_ids[9], node_ids[5]});
        case 4:
            return std::make_shared<Surface6>(
                std::array<ID, 6> {node_ids[2], node_ids[3], node_ids[0], node_ids[9], node_ids[7], node_ids[6]});
        default: return nullptr;    // Invalid surface ID
    }
}
const quadrature::Quadrature& C3D10::integration_scheme() {
    const static quadrature::Quadrature quad {quadrature::DOMAIN_ISO_TET, quadrature::ORDER_QUADRATIC};
    return quad;
}
const quadrature::Quadrature& C3D10::integration_scheme_mass() {
    const static quadrature::Quadrature quad {quadrature::DOMAIN_ISO_TET, quadrature::ORDER_QUARTIC};
    return quad;
}

}    // namespace model
}    // namespace fem