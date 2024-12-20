//
// Created by Luecx on 12.06.2023.
//

#include "c3d8.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface3.h"

namespace fem {
namespace model {

C3D8::C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

StaticMatrix<8, 1> C3D8::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<8, 1>  res {};

    Precision rp = r + 1;
    Precision sp = s + 1;
    Precision tp = t + 1;
    Precision rm = r - 1;
    Precision sm = s - 1;
    Precision tm = t - 1;

    // shape function evaluation. 4 entries for 1) sign, 2) r 3) s 4) t
    for (int n = 0; n < 8; n++) {
        Precision sign = (n == 0 || n == 2 || n == 5 || n == 7) ? -1 : 1;
        Precision r1   = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        Precision s1   = ((n) / 2) % 2 == 1 ? sp : sm;
        Precision t1   = n >= 4 ? tp : tm;

        res(n, 0)      = sign * r1 * s1 * t1;
    }

    res *= 0.125;
    return res;
}

StaticMatrix<8, 3> C3D8::shape_derivative(Precision r, Precision s, Precision t) {
    Precision rp = r + 1;
    Precision sp = s + 1;
    Precision tp = t + 1;
    Precision rm = r - 1;
    Precision sm = s - 1;
    Precision tm = t - 1;

    // shape function evaluation. 4 entries for 1) sign, 2) r 3) s 4) t
    StaticMatrix<8, 4> shape_function {};
    for (int n = 0; n < 8; n++) {
        shape_function(n, 0) = (n == 0 || n == 2 || n == 5 || n == 7) ? -1 : 1;
        shape_function(n, 1) = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        shape_function(n, 2) = ((n) / 2) % 2 == 1 ? sp : sm;
        shape_function(n, 3) = n >= 4 ? tp : tm;
    }

    // containing derivatives of the shape functions h1 -> h8 with respect to r/s/t
    StaticMatrix<8, 3> local_shape_derivative {};
    for (int n = 0; n < 8; n++) {
        local_shape_derivative(n, 0) = shape_function(n, 0) * shape_function(n, 2) * shape_function(n, 3);
        local_shape_derivative(n, 1) = shape_function(n, 0) * shape_function(n, 1) * shape_function(n, 3);
        local_shape_derivative(n, 2) = shape_function(n, 0) * shape_function(n, 1) * shape_function(n, 2);
    }
    local_shape_derivative *= 0.125;
    return local_shape_derivative;
}

StaticMatrix<8, 3> C3D8::node_coords_local() {
    StaticMatrix<8, 3> res {};
    for (int n = 0; n < 8; n++) {
        Precision r1 = ((n + 1) / 2) % 2 == 0 ? -1 : 1;
        Precision s1 = ((n) / 2) % 2 == 0 ? -1 : 1;
        Precision t1 = n >= 4 ? 1 : -1;

        res(n, 0)    = r1;
        res(n, 1)    = s1;
        res(n, 2)    = t1;
    }
    return res;
}
SurfacePtr C3D8::surface(ID surface_id) {
    switch (surface_id) {
        case 1:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[0], node_ids[1], node_ids[2], node_ids[3]});
        case 2:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[4], node_ids[5], node_ids[6], node_ids[7]});
        case 3:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[0], node_ids[1], node_ids[5], node_ids[4]});
        case 4:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[1], node_ids[2], node_ids[6], node_ids[5]});
        case 5:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[2], node_ids[3], node_ids[7], node_ids[6]});
        case 6:
            return std::make_shared<Surface4>(std::array<ID, 4> {node_ids[3], node_ids[0], node_ids[4], node_ids[7]});
        default: return nullptr;    // Invalid surface ID
    }
}
const quadrature::Quadrature& C3D8::integration_scheme() {
    const static quadrature::Quadrature quad {quadrature::DOMAIN_ISO_HEX, quadrature::ORDER_QUADRATIC};
    return quad;
}

}    // namespace model
}    // namespace fem