/**
 * @file surface4.h
 * @brief Defines the Surface4 class for a quadrilateral surface element.
 *
 * @details The Surface4 class represents a 4-node quadrilateral surface element
 *          (S3D4). It inherits from the Surface<4> base class and provides
 *          shape function definitions, shape derivatives, and integration schemes
 *          specific to quadrilateral elements. The class supports operations such
 *          as projection onto the surface, computation of local coordinates, and
 *          surface integration.
 *
 * @date Created on 27.09.2024 by Finn Eggers
 */

#pragma once

#include "surface.h"
#include "../line/line2a.h"
#include "../line/line2b.h"

#include <array>
#include <tuple>

namespace fem::model {

/**
 * @class Surface4
 * @brief Represents a quadrilateral surface element with 4 nodes (S3D4).
 *
 * @details The Surface4 class provides shape function definitions, shape
 *          derivatives, and integration schemes tailored for quadrilateral
 *          elements. It supports operations specific to quadrilaterals,
 *          such as checking element bounds, projecting points onto the
 *          element, and determining the closest point on the element's
 *          boundary.
 */
struct Surface4 : public Surface<4> {

    using Surface<4>::num_nodes;
    using Surface<4>::num_edges;
    using Surface<4>::num_nodes_per_edge;

    /**
     * @brief Constructor for the Surface4 class.
     *
     * @param pNodeIds Array of node IDs corresponding to the surface element nodes.
     */
    Surface4(const std::array<ID, 4>& pNodeIds);

    /**
     * @brief Compute shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<4, 1> Vector of shape function values at (r, s).
     */
    StaticMatrix<4, 1> shape_function(Precision r, Precision s) const override;

    /**
     * @brief Compute shape function derivatives at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<4, 2> Matrix of shape function derivatives at (r, s).
     */
    StaticMatrix<4, 2> shape_derivative(Precision r, Precision s) const override;

    /**
     * @brief Compute second-order derivatives of shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<4, 3> Matrix of second-order shape function derivatives.
     */
    StaticMatrix<4, 3> shape_second_derivative(Precision r, Precision s) const override;

    /**
     * @brief Retrieve the local coordinates of the nodes for the quadrilateral element.
     *
     * @return StaticMatrix<4, 2> Local coordinates of the nodes.
     */
    StaticMatrix<4, 2> node_coords_local() const override;

    /**
     * @brief Get the integration scheme for a quadrilateral element.
     *
     * @return const quadrature::Quadrature& Integration scheme using a linear quadrature rule.
     */
    const quadrature::Quadrature& integration_scheme() const override;

    /**
     * @brief Compute the closest point on the element boundary to a given global point.
     *
     * @param global Global coordinates of the point.
     * @param node_coords Coordinates of the nodes in the global system.
     * @return Vec2 Local coordinates of the closest boundary point.
     */
    Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<4, 3>& node_coords) const override;

    /**
     * @brief Check if a given local coordinate is within the bounds of the element.
     *
     * @param local The local coordinates (r, s) to check.
     * @return bool True if the point is within bounds, false otherwise.
     */
    bool in_bounds(const Vec2& local) const override;
};

}  // namespace fem::model
