/**
 * @file surface_geometry.inl
 * @brief Implements interpolation, geometry and connectivity for finite-element surfaces.
 *
 * These routines provide the topology-independent parts of the isoparametric
 * surface description. They gather nodal coordinates, evaluate mappings and
 * Jacobians, and expose the natural domain used by projection and integration.
 *
 * @see Surface
 * @see SurfaceInterface
 */

#pragma once

namespace fem::model {

/**
 * @brief Constructs a finite-element surface from its global node identifiers.
 *
 * The topology constants are passed to `SurfaceInterface`, while the node
 * identifiers are retained in fixed-size local storage.
 *
 * @param node_ids Global identifiers of the surface nodes.
 */
template<Index N>
Surface<N>::Surface(const std::array<ID, N>& node_ids)
    : SurfaceInterface(num_edges, num_nodes, num_nodes_per_edge),
      nodeIds         (node_ids) {}

/**
 * @brief Evaluates the shape functions using dynamic vector storage.
 *
 * Concrete surfaces implement the fixed-size two-coordinate overload. This
 * adapter exposes those shape functions through the common dynamic interface.
 */
template<Index N>
DynamicVector Surface<N>::shape_function(const Vec2& local) const {
    return DynamicVector(shape_function(local(0), local(1)));
}

/**
 * @brief Evaluates first natural shape-function derivatives dynamically.
 */
template<Index N>
DynamicMatrix Surface<N>::shape_derivative(const Vec2& local) const {
    return DynamicMatrix(shape_derivative(local(0), local(1)));
}

/**
 * @brief Evaluates second natural shape-function derivatives dynamically.
 *
 * Columns contain d2N/dr2, d2N/ds2 and d2N/(dr ds), matching the fixed-size
 * surface implementations.
 */
template<Index N>
DynamicMatrix Surface<N>::shape_second_derivative(const Vec2& local) const {
    return DynamicMatrix(shape_second_derivative(local(0), local(1)));
}

/**
 * @brief Interpolates fixed-size nodal values at natural surface coordinates.
 *
 * The shape-function vector provides the interpolation weights for every
 * surface node. Each column of `nodal_values` is interpolated independently.
 *
 * @tparam M Number of values stored per node.
 *
 * @param nodal_values Nodal values arranged as one node per row.
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Interpolated value vector at `(r,s)`.
 */
template<Index N>
template<int M>
StaticVector<M> Surface<N>::interpolate(
    const StaticMatrix<N, M>& nodal_values,
    Precision                 r,
    Precision                 s
) const {
    return (shape_function(r, s).transpose() * nodal_values).transpose();
}

/**
 * @brief Computes the local-to-global surface Jacobian.
 *
 * Multiplying the global nodal coordinates by the shape-function derivatives
 * produces the two physical tangent vectors `dx/dr` and `dx/ds`.
 */
template<Index N>
StaticMatrix<3, 2> Surface<N>::jacobian(
    const StaticMatrix<N, 3>& node_coords,
    Precision                 r,
    Precision                 s
) const {
    return node_coords.transpose() * shape_derivative(r, s);
}

/**
 * @brief Returns the polygon of the natural element domain.
 */
template<Index N>
SurfaceInterface::Polygon Surface<N>::local_domain_polygon() const {
    if constexpr (N == 3 || N == 6) {
        return Polygon{
            Vec2(Precision(0), Precision(0)),
            Vec2(Precision(1), Precision(0)),
            Vec2(Precision(0), Precision(1))
        };
    } else {
        return Polygon{
            Vec2(Precision(-1), Precision(-1)),
            Vec2(Precision( 1), Precision(-1)),
            Vec2(Precision( 1), Precision( 1)),
            Vec2(Precision(-1), Precision( 1))
        };
    }
}

/**
 * @brief Maps natural surface coordinates into physical space.
 */
template<Index N>
Vec3 Surface<N>::local_to_global(const Vec2& local, const Field& node_coords) const {
    const auto coordinates = node_coords_global(node_coords);
    return interpolate(coordinates, local(0), local(1));
}

/**
 * @brief Returns the natural coordinates of all surface nodes dynamically.
 */
template<Index N>
DynamicMatrix Surface<N>::node_coords_natural() const {
    return node_coords_local();
}

/**
 * @brief Computes the normalized physical surface normal.
 */
template<Index N>
Vec3 Surface<N>::normal(const Field& node_coords, const Vec2& local) const {
    const auto coordinates = node_coords_global(node_coords);
    const auto jac         = jacobian(coordinates, local(0), local(1));

    return jac.col(0).cross(jac.col(1)).normalized();
}

/**
 * @brief Collects the global coordinates of all nodes attached to the surface.
 */
template<Index N>
StaticMatrix<N, 3> Surface<N>::node_coords_global(const Field& node_coords) const {
    StaticMatrix<N, 3> coordinates{};

    for (Index local_id = 0; local_id < N; ++local_id) {
        const Index node_id = static_cast<Index>(nodeIds[local_id]);
        coordinates.row(local_id) = node_coords.row_vec3(node_id).transpose();
    }

    return coordinates;
}

/**
 * @brief Returns a mutable pointer to the first global node identifier.
 */
template<Index N>
ID* Surface<N>::nodes() {
    return nodeIds.data();
}

/**
 * @brief Returns a const pointer to the first global node identifier.
 */
template<Index N>
const ID* Surface<N>::nodes() const {
    return nodeIds.data();
}

} // namespace fem::model
