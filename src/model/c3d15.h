#pragma once

#include "element.h"
#include "element_solid.h"

namespace fem { namespace model {

struct C3D15 : public SolidElement<15>{
    C3D15(ID p_elem_id, const std::array<ID, 15>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_WEDGE, quadrature::ORDER_SUPER_QUADRATIC};
        return quad;
    }

    StaticMatrix<15, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> node_coords_local() override;

};

} }
