#pragma once

#include "element.h"
#include "element_solid.h"
#include "surface/surface8.h"
#include "surface/surface6.h"
#include "surface/surface4.h"
#include "surface/surface3.h"

namespace fem { namespace model {

struct C3D13 : public SolidElement<13>{
    C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_PYRAMID, quadrature::ORDER_QUARTIC};
        return quad;
    }

    SurfacePtr surface(ID surface_id) override {
        return nullptr;
        // raise std::runtime_error("Not implemented");
    }

    StaticMatrix<13, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> node_coords_local() override;

};

} }
