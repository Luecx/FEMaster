/**
 * @file surface4.cpp
 * @brief Implements the four-node bilinear quadrilateral surface element.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#include "surface4.h"

#include "../line/line2a.h"

namespace fem::model {

/**
 * @brief Constructs a four-node quadrilateral surface.
 *
 * The node identifiers follow the natural corner ordering
 * `(-1,-1)`, `(1,-1)`, `(1,1)`, and `(-1,1)`.
 *
 * @param node_ids Global identifiers of the four surface nodes.
 */
Surface4::Surface4(const std::array<ID, 4>& node_ids)
    : Surface<4>(node_ids) {}

/**
 * @brief Evaluates the bilinear quadrilateral shape functions.
 *
 * Each shape function is one at its associated corner and zero at the
 * remaining three corners of the natural element domain.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Shape-function vector evaluated at `(r,s)`.
 */
StaticMatrix<4, 1> Surface4::shape_function(Precision r, Precision s) const {
    StaticMatrix<4, 1> shape;

    // Evaluate the four bilinear corner shape functions following the node
    // Ordering around the natural quadrilateral domain
    shape << Precision(0.25) * (Precision(1) - r) * (Precision(1) - s),
             Precision(0.25) * (Precision(1) + r) * (Precision(1) - s),
             Precision(0.25) * (Precision(1) + r) * (Precision(1) + s),
             Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);

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
StaticMatrix<4, 2> Surface4::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<4, 2> derivative;

    // Each row belongs to one shape function while both columns contain its
    // Derivatives with respect to r and s
    derivative(0, 0) = -Precision(0.25) * (Precision(1) - s);
    derivative(0, 1) = -Precision(0.25) * (Precision(1) - r);

    derivative(1, 0) =  Precision(0.25) * (Precision(1) - s);
    derivative(1, 1) = -Precision(0.25) * (Precision(1) + r);

    derivative(2, 0) =  Precision(0.25) * (Precision(1) + s);
    derivative(2, 1) =  Precision(0.25) * (Precision(1) + r);

    derivative(3, 0) = -Precision(0.25) * (Precision(1) + s);
    derivative(3, 1) =  Precision(0.25) * (Precision(1) - r);

    return derivative;
}

/**
 * @brief Evaluates the second derivatives of the shape functions.
 *
 * Bilinear shape functions have no pure quadratic terms in either natural
 * coordinate. Therefore, `d²N/dr²` and `d²N/ds²` vanish while the mixed
 * derivatives `d²N/(dr ds)` remain constant.
 *
 * The three matrix columns contain `d²N/dr²`, `d²N/ds²`, and
 * `d²N/(dr ds)`, respectively.
 *
 * @param r First natural coordinate; unused because the derivatives are constant.
 * @param s Second natural coordinate; unused because the derivatives are constant.
 *
 * @return Matrix containing the second shape-function derivatives.
 */
StaticMatrix<4, 3> Surface4::shape_second_derivative(Precision r, Precision s) const {
    // The second derivatives are constant, but the coordinates remain part of
    // The common surface interface
    (void) r;
    (void) s;

    StaticMatrix<4, 3> second_derivative;

    // Only the mixed derivatives are non-zero for bilinear interpolation
    second_derivative << Precision(0), Precision(0),  Precision(0.25),
                         Precision(0), Precision(0), -Precision(0.25),
                         Precision(0), Precision(0),  Precision(0.25),
                         Precision(0), Precision(0), -Precision(0.25);

    return second_derivative;
}

/**
 * @brief Returns the node positions in the natural coordinate system.
 *
 * The four nodes are placed at the corners of the square natural domain
 * `[-1,1] x [-1,1]`.
 *
 * @return Matrix containing one local `(r,s)` coordinate per node.
 */
StaticMatrix<4, 2> Surface4::node_coords_local() const {
    StaticMatrix<4, 2> local_coords;

    // The row order has to match the node and shape-function ordering
    local_coords << Precision(-1), Precision(-1),  // node 1
                    Precision( 1), Precision(-1),  // node 2
                    Precision( 1), Precision( 1),  // node 3
                    Precision(-1), Precision( 1);  // node 4

    return local_coords;
}

/**
 * @brief Computes the closest point on the quadrilateral boundary.
 *
 * The four element edges are represented by linear line elements. The global
 * point is projected onto every edge, after which the closest projection is
 * transformed back into the natural coordinates of the quadrilateral.
 *
 * @param global Global point to project onto the element boundary.
 * @param node_coords Global coordinates of the four quadrilateral nodes.
 *
 * @return Natural quadrilateral coordinates of the closest boundary point.
 */
Vec2 Surface4::closest_point_on_boundary(const Vec3&               global,
                                         const StaticMatrix<4, 3>& node_coords) const {
    // Represent all four quadrilateral edges by line elements following the
    // Counter-clockwise surface node ordering
    Line2A edge_01({0, 1});
    Line2A edge_12({1, 2});
    Line2A edge_23({2, 3});
    Line2A edge_30({3, 0});

    // The line-element interface operates on Field storage, so transfer the
    // Supplied fixed-size coordinate matrix into a temporary nodal field
    Field node_field("SURFACE4_BOUNDARY", FieldDomain::NODE, 4, 3);

    for (Index local_id = 0; local_id < 4; ++local_id) {
        for (Dim component = 0; component < 3; ++component) {
            node_field(local_id, component) = node_coords(local_id, component);
        }
    }

    // Project the global point independently onto every quadrilateral edge
    const Precision edge_01_local = edge_01.global_to_local(global, node_field);
    const Precision edge_12_local = edge_12.global_to_local(global, node_field);
    const Precision edge_23_local = edge_23.global_to_local(global, node_field);
    const Precision edge_30_local = edge_30.global_to_local(global, node_field);

    // Map the edge-local projections back into physical coordinates so their
    // Distances to the requested global point can be compared
    const Vec3 point_01 = edge_01.local_to_global(edge_01_local, node_field);
    const Vec3 point_12 = edge_12.local_to_global(edge_12_local, node_field);
    const Vec3 point_23 = edge_23.local_to_global(edge_23_local, node_field);
    const Vec3 point_30 = edge_30.local_to_global(edge_30_local, node_field);

    // Squared distances are sufficient for comparison and avoid unnecessary
    // Square-root evaluations
    const Precision distance_01 = (point_01 - global).squaredNorm();
    const Precision distance_12 = (point_12 - global).squaredNorm();
    const Precision distance_23 = (point_23 - global).squaredNorm();
    const Precision distance_30 = (point_30 - global).squaredNorm();

    // Edge 0-1 follows (r,s) = (p,-1)
    if (distance_01 <= distance_12 && distance_01 <= distance_23 && distance_01 <= distance_30) {
        return {edge_01_local, Precision(-1)};
    }

    // Edge 1-2 follows (r,s) = (1,p)
    if (distance_12 <= distance_01 && distance_12 <= distance_23 && distance_12 <= distance_30) {
        return {Precision(1), edge_12_local};
    }

    // Edge 2-3 runs from r=1 to r=-1 and therefore follows (r,s) = (-p,1)
    if (distance_23 <= distance_01 && distance_23 <= distance_12 && distance_23 <= distance_30) {
        return {-edge_23_local, Precision(1)};
    }

    // Edge 3-0 runs from s=1 to s=-1 and therefore follows (r,s) = (-1,-p)
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
bool Surface4::in_bounds(const Vec2& local) const {
    const Precision r = local(0);
    const Precision s = local(1);

    return r >= Precision(-1) && r <= Precision(1) &&
           s >= Precision(-1) && s <= Precision(1);
}

/**
 * @brief Returns the quadrature rule used for the bilinear quadrilateral.
 *
 * The static quadrature object is constructed only once and reused for all
 * element evaluations.
 *
 * @return Linear-order quadrature rule on the isoparametric square domain.
 */
const math::quadrature::Quadrature& Surface4::integration_scheme() const {
    static const math::quadrature::Quadrature scheme{
        math::quadrature::DOMAIN_ISO_QUAD,
        math::quadrature::ORDER_LINEAR
    };

    return scheme;
}

} // namespace fem::model