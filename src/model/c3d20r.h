#pragma once

#include "element.h"
#include "element_solid.h"
#include "surface/surface8.h"
#include "surface/surface6.h"
#include "surface/surface4.h"
#include "surface/surface3.h"

namespace fem { namespace model {

struct C3D20R : public SolidElement<20>{
    C3D20R(ID p_elem_id, const std::array<ID, 20>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_HEX, quadrature::ORDER_QUADRATIC};
        return quad;
    }


    SurfacePtr surface(ID surface_id) override {
        switch (surface_id) {
            case 1: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 0], node_ids[ 1], node_ids[ 2], node_ids[ 3],
                                                                        node_ids[ 8], node_ids[ 9], node_ids[10], node_ids[11]});
            case 2: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 4], node_ids[ 5], node_ids[ 6], node_ids[ 7],
                                                                        node_ids[12], node_ids[13], node_ids[14], node_ids[15]});
            case 3: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 0], node_ids[ 1], node_ids[ 5], node_ids[ 4],
                                                                        node_ids[ 8], node_ids[17], node_ids[12], node_ids[16]});
            case 4: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 1], node_ids[ 2], node_ids[ 6], node_ids[ 5],
                                                                        node_ids[ 9], node_ids[18], node_ids[13], node_ids[17]});
            case 5: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 2], node_ids[ 3], node_ids[ 7], node_ids[ 6],
                                                                        node_ids[10], node_ids[19], node_ids[14], node_ids[18]});
            case 6: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[ 3], node_ids[ 0], node_ids[ 4], node_ids[ 7],
                                                                        node_ids[11], node_ids[12], node_ids[15], node_ids[19]});
            default: return nullptr;  // Invalid surface ID
        }
    }

    StaticMatrix<20, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> node_coords_local() override;

};

} }
