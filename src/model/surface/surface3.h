/******************************************************************************
 * @file surface3.h
 * @brief Derived class for a triangular surface element with 3 nodes (S3D3).
 *
 * @details This file defines the `Surface3` class, which inherits from the
 *          `SurfaceInterface<3>` base class. It represents a triangular
 *          surface element with three nodes, used primarily for surface
 *          integrals and projecting points onto surfaces. It provides
 *          shape functions, shape function derivatives, and surface
 *          integration schemes tailored for triangular elements.
 *
 *          This class is used in surface element interactions, such as
 *          tie constraints and boundary condition enforcement.
 *
 * @date Created on 27.09.2024 by Finn Eggers
 ******************************************************************************/

#pragma once

#include "../line/line2b.h"
#include "surface.h"

#include <array>
#include <tuple>

namespace fem::model {

/******************************************************************************
 * @class Surface3
 * @brief Represents a triangular surface element with 3 nodes (S3D3).
 *
 * @details This derived class provides shape function definitions, integration
 *          schemes, and boundary checks for a simple triangular surface element.
 *          It supports common surface operations such as projection to the
 *          element, finding local coordinates, and computing surface areas.
 ******************************************************************************/
struct Surface3 : public SurfaceInterface<3> {
    using SurfaceInterface<3>::num_nodes;
    using SurfaceInterface<3>::num_edges;
    using SurfaceInterface<3>::num_nodes_per_edge;

    /**
     * @brief Constructor for the Surface3 class.
     *
     * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
     */
    Surface3(const std::array<ID, 3>& pNodeIds);

    /**
     * @brief Compute shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<3, 1> Vector of shape function values at (r, s).
     */
    StaticMatrix<3, 1> shape_function(Precision r, Precision s) const override;

    /**
     * @brief Compute shape function derivatives at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<3, 2> Matrix of shape function derivatives at (r, s).
     */
    StaticMatrix<3, 2> shape_derivative(Precision r, Precision s) const override;

    /**
     * @brief Compute second-order derivatives of shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<3, 3> Matrix of second-order shape function derivatives.
     */
    StaticMatrix<3, 3> shape_second_derivative(Precision r, Precision s) const override;

    /**
     * @brief Retrieve the local coordinates of the nodes for the triangular element.
     *
     * @return StaticMatrix<3, 2> Local coordinates of the nodes.
     */
    StaticMatrix<3, 2> node_coords_local() const override;

    /**
     * @brief Get the integration scheme for a triangular element.
     *
     * @return const quadrature::Quadrature& Integration scheme using a linear quadrature rule.
     */
    const quadrature::Quadrature& integration_scheme() const override;

    /**
     * @brief Compute the closest point on the element boundary to a given point.
     *
     * @param global Global coordinates of the point.
     * @param node_coords Coordinates of the nodes in the global system.
     * @return Vec2 Local coordinates of the closest boundary point.
     */
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<3, 3>& node_coords) const override;

    /**
     * @brief Check if a given local coordinate is within the bounds of the element.
     *
     * @param local The local coordinates (r, s) to check.
     * @return bool True if the point is within bounds, false otherwise.
     */
    bool in_bounds(const Vec2& local) const override;
};

}    // namespace fem::model
