/**
 * @file surface6.cpp
 * @brief Implements the six-node quadratic triangular surface element.
 *
 * @author Finn Eggers
 * @date 01.10.2024
 */

#include "surface6.h"

#include "../line/line3b.h"

namespace fem::model {

/**
 * @brief Constructs a six-node quadratic triangular surface.
 *
 * The first three node identifiers correspond to the triangle corners. The
 * remaining identifiers correspond to the midside nodes on edges `0-1`,
 * `1-2`, and `2-0`.
 *
 * @param node_ids Global identifiers of the six surface nodes.
 */
Surface6::Surface6(const std::array<ID, 6>& node_ids)
    : Surface<6>(node_ids) {}

/**
 * @brief Evaluates the quadratic triangular shape functions.
 *
 * The interpolation uses the barycentric coordinates
 * `L1 = 1 - r - s`, `L2 = r`, and `L3 = s`. The first three shape
 * functions belong to the corner nodes, while the remaining functions belong
 * to the midside nodes.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Shape-function vector evaluated at `(r,s)`.
 */
StaticMatrix<6, 1> Surface6::shape_function(Precision r, Precision s) const {
    // Define the first barycentric coordinate once because it occurs in
    // Several corner and midside shape functions
    const Precision t = Precision(1) - r - s;

    StaticMatrix<6, 1> shape;

    // The first three entries interpolate the corner nodes while the final
    // Three entries interpolate the midside nodes
    shape << t * (Precision(2) * t - Precision(1)),
             r * (Precision(2) * r - Precision(1)),
             s * (Precision(2) * s - Precision(1)),
             Precision(4) * r * t,
             Precision(4) * r * s,
             Precision(4) * s * t;

    return shape;
}

/**
 * @brief Evaluates the first derivatives of the shape functions.
 *
 * The first and second matrix columns contain `dN/dr` and `dN/ds`,
 * respectively.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Matrix containing the first shape-function derivatives.
 */
StaticMatrix<6, 2> Surface6::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<6, 2> derivative;

    // Each row belongs to one quadratic shape function while the columns
    // Contain differentiation with respect to r and s
    derivative(0, 0) = -Precision(3) + Precision(4) * (r + s);
    derivative(0, 1) = -Precision(3) + Precision(4) * (r + s);

    derivative(1, 0) = Precision(4) * r - Precision(1);
    derivative(1, 1) = Precision(0);

    derivative(2, 0) = Precision(0);
    derivative(2, 1) = Precision(4) * s - Precision(1);

    derivative(3, 0) = Precision(4) - Precision(8) * r - Precision(4) * s;
    derivative(3, 1) = -Precision(4) * r;

    derivative(4, 0) = Precision(4) * s;
    derivative(4, 1) = Precision(4) * r;

    derivative(5, 0) = -Precision(4) * s;
    derivative(5, 1) = Precision(4) - Precision(4) * r - Precision(8) * s;

    return derivative;
}

/**
 * @brief Evaluates the second derivatives of the shape functions.
 *
 * Because all shape functions are quadratic, their second derivatives are
 * constant throughout the natural element domain.
 *
 * The three matrix columns contain `d²N/dr²`, `d²N/ds²`, and
 * `d²N/(dr ds)`, respectively.
 *
 * @param r First natural coordinate; unused because the derivatives are constant.
 * @param s Second natural coordinate; unused because the derivatives are constant.
 *
 * @return Matrix containing the second shape-function derivatives.
 */
StaticMatrix<6, 3> Surface6::shape_second_derivative(Precision r, Precision s) const {
    // The second derivatives are constant, but the coordinates remain part of
    // The common surface interface
    (void) r;
    (void) s;

    StaticMatrix<6, 3> second_derivative;

    // Each row contains the two pure and the mixed second derivative of one
    // Quadratic triangular shape function
    second_derivative <<  Precision(4),  Precision(4),  Precision(4),
                          Precision(4),  Precision(0),  Precision(0),
                          Precision(0),  Precision(4),  Precision(0),
                         -Precision(8),  Precision(0), -Precision(4),
                          Precision(0),  Precision(0),  Precision(4),
                          Precision(0), -Precision(8), -Precision(4);

    return second_derivative;
}

/**
 * @brief Returns the node positions in the natural coordinate system.
 *
 * The first three rows contain the triangle corners. The remaining rows
 * contain the midside positions of edges `0-1`, `1-2`, and `2-0`.
 *
 * @return Matrix containing one local `(r,s)` coordinate per node.
 */
StaticMatrix<6, 2> Surface6::node_coords_local() const {
    StaticMatrix<6, 2> local_coords;

    // The row order has to match the node and shape-function ordering
    local_coords << Precision(0.0), Precision(0.0),  // corner node 1
                    Precision(1.0), Precision(0.0),  // corner node 2
                    Precision(0.0), Precision(1.0),  // corner node 3
                    Precision(0.5), Precision(0.0),  // midside node 4
                    Precision(0.5), Precision(0.5),  // midside node 5
                    Precision(0.0), Precision(0.5);  // midside node 6

    return local_coords;
}

/**
 * @brief Computes the closest point on the quadratic triangular boundary.
 *
 * The three curved triangle edges are represented by quadratic line elements.
 * The global point is projected onto every edge, after which the closest
 * projection is transformed back into the natural coordinates of the triangle.
 *
 * @param global Global point to project onto the element boundary.
 * @param node_coords Global coordinates of the six triangle nodes.
 *
 * @return Natural triangle coordinates of the closest boundary point.
 */
Vec2 Surface6::closest_point_on_boundary(const Vec3&               global,
                                         const StaticMatrix<6, 3>& node_coords) const {
    // Represent the three quadratic triangle edges using both corner nodes and
    // The corresponding midside node
    Line3B edge_01({0, 1, 3});
    Line3B edge_12({1, 2, 4});
    Line3B edge_20({2, 0, 5});

    // The line-element interface operates on Field storage, so transfer the
    // Supplied fixed-size coordinate matrix into a temporary nodal field
    Field node_field("SURFACE6_BOUNDARY", FieldDomain::NODE, 6, 3);

    for (Index local_id = 0; local_id < 6; ++local_id) {
        for (Dim component = 0; component < 3; ++component) {
            node_field(local_id, component) = node_coords(local_id, component);
        }
    }

    // Project the global point independently onto every curved triangle edge
    const Precision edge_01_local = edge_01.global_to_local(global, node_field);
    const Precision edge_12_local = edge_12.global_to_local(global, node_field);
    const Precision edge_20_local = edge_20.global_to_local(global, node_field);

    // Map the edge-local projections back into physical coordinates so their
    // Distances to the requested global point can be compared
    const Vec3 point_01 = edge_01.local_to_global(edge_01_local, node_field);
    const Vec3 point_12 = edge_12.local_to_global(edge_12_local, node_field);
    const Vec3 point_20 = edge_20.local_to_global(edge_20_local, node_field);

    // Squared distances are sufficient for comparison and avoid unnecessary
    // Square-root evaluations
    const Precision distance_01 = (point_01 - global).squaredNorm();
    const Precision distance_12 = (point_12 - global).squaredNorm();
    const Precision distance_20 = (point_20 - global).squaredNorm();

    // Edge 0-1 follows (r,s) = (p,0)
    if (distance_01 <= distance_12 && distance_01 <= distance_20) {
        return {edge_01_local, Precision(0)};
    }

    // Edge 2-0 follows (r,s) = (0,1-p) because its line orientation starts at
    // Node 2 and ends at node 0
    if (distance_20 <= distance_01 && distance_20 <= distance_12) {
        return {Precision(0), Precision(1) - edge_20_local};
    }

    // Edge 1-2 follows (r,s) = (1-p,p)
    return {Precision(1) - edge_12_local, edge_12_local};
}

/**
 * @brief Checks whether natural coordinates lie inside the reference triangle.
 *
 * A point belongs to the triangular domain when both coordinates are
 * non-negative and their sum does not exceed one.
 *
 * @param local Natural coordinates `(r,s)` to test.
 *
 * @return `true` when the point lies inside or on the triangle boundary.
 */
bool Surface6::in_bounds(const Vec2& local) const {
    const Precision r = local(0);
    const Precision s = local(1);

    return r >= Precision(0) &&
           s >= Precision(0) &&
           r + s <= Precision(1);
}

/**
 * @brief Returns the quadrature rule used for the quadratic triangle.
 *
 * The quadratic rule accounts for the higher interpolation order of the
 * six-node surface. The static object is constructed only once and reused for
 * all evaluations.
 *
 * @return Quadratic-order quadrature rule on the isoparametric triangle.
 */
const math::quadrature::Quadrature& Surface6::integration_scheme() const {
    static const math::quadrature::Quadrature scheme{
        math::quadrature::DOMAIN_ISO_TRI,
        math::quadrature::ORDER_QUADRATIC
    };

    return scheme;
}

} // namespace fem::model