/******************************************************************************
 * @file C3D4.h
 * @brief C3D4.h defines a 4-node tetrahedral element (C3D4) for solid
 * mechanics within a finite element model (FEM). This element computes 
 * shape functions, their derivatives, and uses a suitable quadrature scheme 
 * for integration.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de> 
 * all rights reserved
 * @date Created on 27.08.2024
 *
 ******************************************************************************/

#pragma once

#include "element_solid.h"
#include "../geometry/surface/surface8.h"
#include "../geometry/surface/surface6.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface3.h"

namespace fem { namespace model {

/******************************************************************************
 * C3D4 class
 * This class defines a 4-node tetrahedral solid element used in finite 
 * element analysis. It provides methods to compute shape functions, their 
 * derivatives, and handle integration with a linear tetrahedral quadrature 
 * scheme.
 ******************************************************************************/
struct C3D4 : public SolidElement<4> {

    //-------------------------------------------------------------------------
    // Constructor
    //-------------------------------------------------------------------------
    /**
     * @brief Constructs a C3D4 element with a given element ID and node IDs.
     *
     * @param pElemId The unique ID of the element.
     * @param pNodeIds Array containing IDs of the 4 nodes.
     */
    C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds);

    //-------------------------------------------------------------------------
    // Shape Function
    //-------------------------------------------------------------------------
    /**
     * @brief Computes the shape functions of the C3D4 element at the given
     * local coordinates (r, s, t).
     *
     * @param r The local coordinate in the r-direction.
     * @param s The local coordinate in the s-direction.
     * @param t The local coordinate in the t-direction.
     * @return StaticMatrix<4, 1> The evaluated shape function values.
     */
    StaticMatrix<4, 1> shape_function(Precision r, Precision s, Precision t) override;

    //-------------------------------------------------------------------------
    // Shape Function Derivatives
    //-------------------------------------------------------------------------
    /**
     * @brief Computes the derivatives of the shape functions with respect to
     * the local coordinates (r, s, t).
     *
     * @param r The local coordinate in the r-direction.
     * @param s The local coordinate in the s-direction.
     * @param t The local coordinate in the t-direction.
     * @return StaticMatrix<4, 3> The derivatives of the shape functions with 
     * respect to (r, s, t).
     */
    StaticMatrix<4, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    //-------------------------------------------------------------------------
    // Node Local Coordinates
    //-------------------------------------------------------------------------
    /**
     * @brief Returns the local coordinates of the nodes of the C3D4 element.
     *
     * @return StaticMatrix<4, 3> The local coordinates of the element's nodes.
     */
    StaticMatrix<4, 3> node_coords_local() override;

    //-------------------------------------------------------------------------
    // Integration Scheme
    //-------------------------------------------------------------------------
    /**
     * @brief Returns the quadrature integration scheme for the C3D4 element,
     * which uses a linear tetrahedral quadrature rule.
     *
     * @return const quadrature::Quadrature& The quadrature rule to be used
     * for integration over the element's domain.
     */
    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_TET, quadrature::ORDER_LINEAR};
        return quad;
    }

    SurfacePtr surface(ID surface_id) override {
        switch (surface_id) {
            case 1: return std::make_unique<Surface3>(std::array<ID, 3>{node_ids[ 0], node_ids[ 1], node_ids[ 2]});
            case 2: return std::make_unique<Surface3>(std::array<ID, 3>{node_ids[ 0], node_ids[ 3], node_ids[ 1]});
            case 3: return std::make_unique<Surface3>(std::array<ID, 3>{node_ids[ 1], node_ids[ 3], node_ids[ 2]});
            case 4: return std::make_unique<Surface3>(std::array<ID, 3>{node_ids[ 2], node_ids[ 3], node_ids[ 0]});
            default: return nullptr;  // Invalid surface ID
        }
    }
};

} } // namespace fem::model
