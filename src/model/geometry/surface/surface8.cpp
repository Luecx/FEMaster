/**
 * @file surface8.cpp
 * @brief Implements the eight-node quadratic quadrilateral surface element.
 *
 * @author Finn Eggers
 * @date 01.10.2024
 */

#include "surface8.h"

#include "../line/line3a.h"

namespace fem::model {

/**
 * @brief Constructs an eight-node quadratic quadrilateral surface.
 *
 * The first four node identifiers correspond to the quadrilateral corners.
 * The remaining identifiers correspond to the midside nodes on edges `0-1`,
 * `1-2`, `2-3`, and `3-0`.
 *
 * @param node_ids Global identifiers of the eight surface nodes.
 */
Surface8::Surface8(const std::array<ID, 8>& node_ids)
    : Surface<8>(node_ids) {}

/**
 * @brief Evaluates the eight-node serendipity shape functions.
 *
 * The first four functions interpolate the corner nodes, while the remaining
 * four functions interpolate the midside nodes of the natural square domain.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Shape-function vector evaluated at `(r,s)`.
 */
StaticMatrix<8, 1> Surface8::shape_function(Precision r, Precision s) const {
    StaticMatrix<8, 1> shape;

    // corner shape functions
    shape(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s) * (-Precision(1) - r - s);
    shape(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s) * (-Precision(1) + r - s);
    shape(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s) * (-Precision(1) + r + s);
    shape(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s) * (-Precision(1) - r + s);

    // midside shape functions
    shape(4) = Precision(0.5) * (Precision(1) - r * r) * (Precision(1) - s);
    shape(5) = Precision(0.5) * (Precision(1) + r) * (Precision(1) - s * s);
    shape(6) = Precision(0.5) * (Precision(1) - r * r) * (Precision(1) + s);
    shape(7) = Precision(0.5) * (Precision(1) - r) * (Precision(1) - s * s);

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
StaticMatrix<8, 2> Surface8::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<8, 2> derivative;

    // derivatives of the four corner shape functions
    derivative(0, 0) = Precision(0.25) * (-Precision(2) * r - s) * (s - Precision(1));
    derivative(0, 1) = Precision(0.25) * (-r - Precision(2) * s) * (r - Precision(1));

    derivative(1, 0) = Precision(0.25) * (-Precision(2) * r + s) * (s - Precision(1));
    derivative(1, 1) = Precision(0.25) * (-r + Precision(2) * s) * (r + Precision(1));

    derivative(2, 0) = Precision(0.25) * (Precision(2) * r + s) * (s + Precision(1));
    derivative(2, 1) = Precision(0.25) * (r + Precision(1)) * (r + Precision(2) * s);

    derivative(3, 0) = Precision(0.25) * (Precision(2) * r - s) * (s + Precision(1));
    derivative(3, 1) = Precision(0.25) * (r - Precision(1)) * (r - Precision(2) * s);

    // derivatives of the four midside shape functions
    derivative(4, 0) = r * (s - Precision(1));
    derivative(4, 1) = Precision(0.5) * (r * r - Precision(1));

    derivative(5, 0) = Precision(0.5) * (Precision(1) - s * s);
    derivative(5, 1) = -s * (r + Precision(1));

    derivative(6, 0) = -r * (Precision(1) + s);
    derivative(6, 1) = Precision(0.5) * (Precision(1) - r * r);

    derivative(7, 0) = Precision(0.5) * (s * s - Precision(1));
    derivative(7, 1) = s * (r - Precision(1));

    return derivative;
}

/**
 * @brief Evaluates the second derivatives of the shape functions.
 *
 * The three matrix columns contain `d²N/dr²`, `d²N/ds²`, and
 * `d²N/(dr ds)`, respectively.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Matrix containing the second shape-function derivatives.
 */
StaticMatrix<8, 3> Surface8::shape_second_derivative(Precision r, Precision s) const {
    StaticMatrix<8, 3> second_derivative;

    // second derivatives of the four corner shape functions
    second_derivative(0, 0) = Precision(0.5) - Precision(0.5) * s;
    second_derivative(0, 1) = Precision(0.5) - Precision(0.5) * r;
    second_derivative(0, 2) = -Precision(0.5) * r - Precision(0.5) * s + Precision(0.25);

    second_derivative(1, 0) = Precision(0.5) - Precision(0.5) * s;
    second_derivative(1, 1) = Precision(0.5) * r + Precision(0.5);
    second_derivative(1, 2) = -Precision(0.5) * r + Precision(0.5) * s - Precision(0.25);

    second_derivative(2, 0) = Precision(0.5) * s + Precision(0.5);
    second_derivative(2, 1) = Precision(0.5) * r + Precision(0.5);
    second_derivative(2, 2) = Precision(0.5) * r + Precision(0.5) * s + Precision(0.25);

    second_derivative(3, 0) = Precision(0.5) * s + Precision(0.5);
    second_derivative(3, 1) = Precision(0.5) - Precision(0.5) * r;
    second_derivative(3, 2) = Precision(0.5) * r - Precision(0.5) * s - Precision(0.25);

    // second derivatives of the four midside shape functions
    second_derivative(4, 0) = s - Precision(1);
    second_derivative(4, 1) = Precision(0);
    second_derivative(4, 2) = r;

    second_derivative(5, 0) = Precision(0);
    second_derivative(5, 1) = -r - Precision(1);
    second_derivative(5, 2) = -s;

    second_derivative(6, 0) = -s - Precision(1);
    second_derivative(6, 1) = Precision(0);
    second_derivative(6, 2) = -r;

    second_derivative(7, 0) = Precision(0);
    second_derivative(7, 1) = r - Precision(1);
    second_derivative(7, 2) = s;

    return second_derivative;
}

/**
 * @brief Returns the node positions in the natural coordinate system.
 *
 * The first four rows contain the corner coordinates. The remaining rows
 * contain the midside coordinates of edges `0-1`, `1-2`, `2-3`, and `3-0`.
 *
 * @return Matrix containing one local `(r,s)` coordinate per node.
 */
StaticMatrix<8, 2> Surface8::node_coords_local() const {
    StaticMatrix<8, 2> local_coords;

    // the row order has to match the node and shape-function ordering
    local_coords << Precision(-1), Precision(-1),  // corner node 1
                    Precision( 1), Precision(-1),  // corner node 2
                    Precision( 1), Precision( 1),  // corner node 3
                    Precision(-1), Precision( 1),  // corner node 4
                    Precision( 0), Precision(-1),  // midside node 5
                    Precision( 1), Precision( 0),  // midside node 6
                    Precision( 0), Precision( 1),  // midside node 7
                    Precision(-1), Precision( 0);  // midside node 8

    return local_coords;
}

/**
 * @brief Computes the closest point on the quadratic quadrilateral boundary.
 *
 * The four curved edges are represented by quadratic line elements. The
 * global point is projected onto every edge, after which the closest projection
 * is transformed back into the natural coordinates of the quadrilateral.
 *
 * @param global Global point to project onto the element boundary.
 * @param node_coords Global coordinates of the eight quadrilateral nodes.
 *
 * @return Natural quadrilateral coordinates of the closest boundary point.
 */
Vec2 Surface8::closest_point_on_boundary(const Vec3&               global,
                                         const StaticMatrix<8, 3>& node_coords) const {
    // represent the four quadratic quadrilateral edges using both corner nodes
    // and the corresponding midside node
    Line3A edge_01({0, 1, 4});
    Line3A edge_12({1, 2, 5});
    Line3A edge_23({2, 3, 6});
    Line3A edge_30({3, 0, 7});

    // the line-element interface operates on Field storage, so transfer the
    // supplied fixed-size coordinate matrix into a temporary nodal field
    Field node_field("SURFACE8_BOUNDARY", FieldDomain::NODE, 8, 3);

    for (Index local_id = 0; local_id < 8; ++local_id) {
        for (Dim component = 0; component < 3; ++component) {
            node_field(local_id, component) = node_coords(local_id, component);
        }
    }

    // project the global point independently onto every curved edge
    const Precision edge_01_local = edge_01.global_to_local(global, node_field);
    const Precision edge_12_local = edge_12.global_to_local(global, node_field);
    const Precision edge_23_local = edge_23.global_to_local(global, node_field);
    const Precision edge_30_local = edge_30.global_to_local(global, node_field);

    // map the edge-local projections back into physical coordinates so their
    // distances to the requested global point can be compared
    const Vec3 point_01 = edge_01.local_to_global(edge_01_local, node_field);
    const Vec3 point_12 = edge_12.local_to_global(edge_12_local, node_field);
    const Vec3 point_23 = edge_23.local_to_global(edge_23_local, node_field);
    const Vec3 point_30 = edge_30.local_to_global(edge_30_local, node_field);

    // squared distances are sufficient for comparison and avoid unnecessary
    // square-root evaluations
    const Precision distance_01 = (point_01 - global).squaredNorm();
    const Precision distance_12 = (point_12 - global).squaredNorm();
    const Precision distance_23 = (point_23 - global).squaredNorm();
    const Precision distance_30 = (point_30 - global).squaredNorm();

    // edge 0-1 follows (r,s) = (p,-1)
    if (distance_01 <= distance_12 && distance_01 <= distance_23 && distance_01 <= distance_30) {
        return {edge_01_local, Precision(-1)};
    }

    // edge 1-2 follows (r,s) = (1,p)
    if (distance_12 <= distance_01 && distance_12 <= distance_23 && distance_12 <= distance_30) {
        return {Precision(1), edge_12_local};
    }

    // edge 2-3 runs from r=1 to r=-1 and therefore follows (r,s) = (-p,1)
    if (distance_23 <= distance_01 && distance_23 <= distance_12 && distance_23 <= distance_30) {
        return {-edge_23_local, Precision(1)};
    }

    // edge 3-0 runs from s=1 to s=-1 and therefore follows (r,s) = (-1,-p)
    return {Precision(-1), -edge_30_local};
}

/**
 * @brief Checks whether natural coordinates lie inside the reference square.
 *
 * A point belongs to the quadrilateral domain when both natural coordinates
 * lie in the closed interval `[-1,1]`.
 *
 * @param local Natural coordinates `(r,s)` to test.
 *
 * @return `true` when the point lies inside or on the element boundary.
 */
bool Surface8::in_bounds(const Vec2& local) const {
    const Precision r = local(0);
    const Precision s = local(1);

    return r >= Precision(-1) && r <= Precision(1) &&
           s >= Precision(-1) && s <= Precision(1);
}

/**
 * @brief Returns the quadrature rule used for the quadratic quadrilateral.
 *
 * The quadratic rule accounts for the higher interpolation order of the
 * eight-node surface. The static object is constructed only once and reused
 * for all evaluations.
 *
 * @return Quadratic-order quadrature rule on the isoparametric square domain.
 */
const math::quadrature::Quadrature& Surface8::integration_scheme() const {
    static const math::quadrature::Quadrature scheme{
        math::quadrature::DOMAIN_ISO_QUAD,
        math::quadrature::ORDER_QUADRATIC
    };

    return scheme;
}

} // namespace fem::model