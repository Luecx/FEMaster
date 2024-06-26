#pragma once

#include "element.h"
#include "element_solid.h"

namespace fem { namespace model {

struct C3D8 : public SolidElement<8>{

    C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds);

    StaticMatrix<8, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<8, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<8, 3> node_coords_local() override;

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_HEX, quadrature::ORDER_QUADRATIC};
        return quad;
    }

};

} }
