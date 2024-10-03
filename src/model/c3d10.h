#pragma once

#include "element.h"
#include "element_solid.h"
#include "surface/surface8.h"
#include "surface/surface6.h"
#include "surface/surface4.h"
#include "surface/surface3.h"

namespace fem { namespace model {

struct C3D10 : public SolidElement<10>{

    C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds);

    StaticMatrix<10, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<10, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<10, 3> node_coords_local() override;

    SurfacePtr surface(ID surface_id) {
        switch (surface_id) {
            case 1: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[0], node_ids[1], node_ids[2], node_ids[4], node_ids[5], node_ids[6]});
            case 2: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[0], node_ids[3], node_ids[1], node_ids[7], node_ids[8], node_ids[4]});
            case 3: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[1], node_ids[3], node_ids[2], node_ids[8], node_ids[9], node_ids[5]});
            case 4: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[2], node_ids[3], node_ids[0], node_ids[9], node_ids[7], node_ids[6]});
            default: return nullptr;  // Invalid surface ID
        }
    }

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_TET, quadrature::ORDER_QUADRATIC};
        return quad;
    }
};

} }
