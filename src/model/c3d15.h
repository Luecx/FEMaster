#pragma once

#include "element.h"
#include "element_solid.h"
#include "surface/surface8.h"
#include "surface/surface6.h"
#include "surface/surface4.h"
#include "surface/surface3.h"

namespace fem { namespace model {

struct C3D15 : public SolidElement<15>{
    C3D15(ID p_elem_id, const std::array<ID, 15>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_WEDGE, quadrature::ORDER_SUPER_QUADRATIC};
        return quad;
    }

    SurfacePtr surface(ID surface_id) override {
        // C3D6 (6-node wedge element): Triangular faces have 3 nodes, quadrilateral faces have 4 nodes
        switch (surface_id) {
            case 1: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[0], node_ids[1], node_ids[2],
                                                                        node_ids[6], node_ids[7], node_ids[8]});  // Face 1: Triangle
            case 2: return std::make_shared<Surface6>(std::array<ID, 6>{node_ids[3], node_ids[4], node_ids[5],
                                                                         node_ids[10], node_ids[11], node_ids[12]});  // Face 2: Triangle
            case 3: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[0], node_ids[1], node_ids[4], node_ids[3],
                                                                         node_ids[6], node_ids[13], node_ids[9], node_ids[12]});  // Face 3: Quadrilateral
            case 4: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[1], node_ids[2], node_ids[5], node_ids[4],
                                                                         node_ids[7], node_ids[14], node_ids[10], node_ids[13]});  // Face 4: Quadrilateral
            case 5: return std::make_shared<Surface8>(std::array<ID, 8>{node_ids[2], node_ids[0], node_ids[3], node_ids[5],
                                                                         node_ids[8], node_ids[12], node_ids[11], node_ids[14]});  // Face 5: Quadrilateral
            default: return nullptr;  // Invalid surface ID
        }
    }

    StaticMatrix<15, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> node_coords_local() override;

};

} }
