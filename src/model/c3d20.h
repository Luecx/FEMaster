#pragma once

#include "element.h"
#include "element_solid.h"

namespace fem { namespace model {

struct C3D20 : public SolidElement<20>{
    C3D20(ID p_elem_id, const std::array<ID, 20>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_HEX, quadrature::ORDER_QUARTIC};
        return quad;
    }

    StaticMatrix<20, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> node_coords_local() override;

};

} }
