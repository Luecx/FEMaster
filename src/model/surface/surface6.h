/******************************************************************************
 * @file surface6.h
 * @brief Defines the Surface6 class for a 6-node triangular surface element.
 *
 * @details This class provides the implementation of a quadratic triangular
 *          surface element used in finite element analysis. It implements
 *          shape functions, shape function derivatives, and local coordinates
 *          mapping specific to 6-node triangles. This element type is often
 *          used for higher-order surface approximations in FEM.
 *
 * @date Created on 01.10.2024 by Finn Eggers
 ******************************************************************************/

#pragma once

#include "surface.h"

namespace fem::model {

/******************************************************************************
 * @class Surface6
 * @brief Represents a triangular surface element with 6 nodes (S3D6).
 *
 * @details The Surface6 class is derived from the `SurfaceInterface<6>` base
 *          class. It implements shape functions, derivatives, and boundary
 *          checks for a quadratic triangular surface element. The class
 *          supports operations such as projecting points onto the surface,
 *          computing local coordinates, and performing surface integration
 *          using a higher-order quadrature scheme.
 ******************************************************************************/
struct Surface6 : public SurfaceInterface<6> {
    using SurfaceInterface<6>::num_nodes;
    using SurfaceInterface<6>::num_edges;
    using SurfaceInterface<6>::num_nodes_per_edge;

    /**
     * @brief Constructor for the Surface6 class.
     *
     * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
     */
    Surface6(const std::array<ID, 6>& pNodeIds);

    /**
     * @brief Compute shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<6, 1> Vector of shape function values at (r, s).
     */
    StaticMatrix<6, 1> shape_function(Precision r, Precision s) const override;

    /**
     * @brief Compute shape function derivatives at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<6, 2> Matrix of shape function derivatives at (r, s).
     */
    StaticMatrix<6, 2> shape_derivative(Precision r, Precision s) const override;

    /**
     * @brief Compute second-order derivatives of shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<6, 3> Matrix of second-order shape function derivatives.
     */
    StaticMatrix<6, 3> shape_second_derivative(Precision r, Precision s) const override;

    /**
     * @brief Retrieve the local coordinates of the nodes for the 6-node triangular element.
     *
     * @return StaticMatrix<6, 2> Local coordinates of the nodes.
     */
    StaticMatrix<6, 2> node_coords_local() const override;

    /**
     * @brief Get the integration scheme for a 6-node triangular element.
     *
     * @return const quadrature::Quadrature& Integration scheme using a quadratic quadrature rule.
     */
    const quadrature::Quadrature& integration_scheme() const override;

    /**
     * @brief Compute the closest point on the element boundary to a given global point.
     *
     * @param global Global coordinates of the point.
     * @param node_coords Coordinates of the nodes in the global system.
     * @return Vec2 Local coordinates of the closest boundary point.
     */
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<6, 3>& node_coords) const override;

    /**
     * @brief Check if a given local coordinate is within the bounds of the element.
     *
     * @param local The local coordinates (r, s) to check.
     * @return bool True if the point is within bounds, false otherwise.
     */
    bool in_bounds(const Vec2& local) const override;
};

}  // namespace fem::model
