#pragma once

#include "element.h"
#include "element_solid.h"
#include "surface/surface8.h"
#include "surface/surface6.h"
#include "surface/surface4.h"
#include "surface/surface3.h"

namespace fem { namespace model {

struct C3D8 : public SolidElement<8>{

    C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds);

    StaticMatrix<8, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<8, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<8, 3> node_coords_local() override;

    SurfacePtr surface(ID surface_id) {
        switch (surface_id) {
            case 1: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 0], node_ids[ 1], node_ids[ 2], node_ids[ 3]});
            case 2: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 4], node_ids[ 5], node_ids[ 6], node_ids[ 7]});
            case 3: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 0], node_ids[ 1], node_ids[ 5], node_ids[ 4]});
            case 4: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 1], node_ids[ 2], node_ids[ 6], node_ids[ 5]});
            case 5: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 2], node_ids[ 3], node_ids[ 7], node_ids[ 6]});
            case 6: return std::make_shared<Surface4>(std::array<ID, 4>{node_ids[ 3], node_ids[ 0], node_ids[ 4], node_ids[ 7]});
            default: return nullptr;  // Invalid surface ID
        }
    }

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_HEX, quadrature::ORDER_QUADRATIC};
        return quad;
    }

};

} }
