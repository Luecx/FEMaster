/**
 * @file surface3.cpp
 * @brief Implements the three-node linear triangular surface element.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#include "surface3.h"

#include "../line/line2b.h"

namespace fem::model {

/**
 * @brief Constructs a three-node triangular surface.
 *
 * The node identifiers follow the local ordering
 * `(0,0)`, `(1,0)`, and `(0,1)`.
 *
 * @param node_ids Global identifiers of the three surface nodes.
 */
Surface3::Surface3(const std::array<ID, 3>& node_ids)
    : Surface<3>(node_ids) {}

/**
 * @brief Evaluates the linear triangular shape functions.
 *
 * The natural coordinates correspond to the barycentric coordinates
 * `L2 = r`, `L3 = s`, and `L1 = 1 - r - s`.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 *
 * @return Shape-function vector evaluated at `(r,s)`.
 */
StaticMatrix<3, 1> Surface3::shape_function(Precision r, Precision s) const {
    StaticMatrix<3, 1> shape;

    // the linear shape functions are identical to the barycentric coordinates
    // of the three-node reference triangle
    shape << Precision(1) - r - s,
             r,
             s;

    return shape;
}

/**
 * @brief Evaluates the first derivatives of the shape functions.
 *
 * Because all shape functions are linear in `r` and `s`, their first
 * derivatives are constant throughout the element.
 *
 * The derivative columns contain `dN/dr` and `dN/ds`, respectively.
 *
 * @param r First natural coordinate; unused for the linear element.
 * @param s Second natural coordinate; unused for the linear element.
 *
 * @return Matrix containing the first shape-function derivatives.
 */
StaticMatrix<3, 2> Surface3::shape_derivative(Precision r, Precision s) const {
    // the derivatives are constant, but the coordinates remain part of the
    // common surface interface
    (void) r;
    (void) s;

    StaticMatrix<3, 2> derivative;

    // each row belongs to one shape function while the columns contain
    // differentiation with respect to r and s
    derivative << Precision(-1), Precision(-1),
                   Precision( 1), Precision( 0),
                   Precision( 0), Precision( 1);

    return derivative;
}

/**
 * @brief Evaluates the second derivatives of the shape functions.
 *
 * Linear triangular shape functions have no quadratic terms. Consequently,
 * all derivatives `d²N/dr²`, `d²N/ds²`, and `d²N/(dr ds)` vanish.
 *
 * @param r First natural coordinate; unused for the linear element.
 * @param s Second natural coordinate; unused for the linear element.
 *
 * @return Zero matrix containing the second shape-function derivatives.
 */
StaticMatrix<3, 3> Surface3::shape_second_derivative(Precision r, Precision s) const {
    // the second derivatives are constant zeros, but the coordinates remain
    // part of the common surface interface
    (void) r;
    (void) s;

    return StaticMatrix<3, 3>::Zero();
}

/**
 * @brief Returns the node positions in the natural coordinate system.
 *
 * The three corner nodes define the reference triangle with the local domain
 * `r >= 0`, `s >= 0`, and `r + s <= 1`.
 *
 * @return Matrix containing one local `(r,s)` coordinate per node.
 */
StaticMatrix<3, 2> Surface3::node_coords_local() const {
    StaticMatrix<3, 2> local_coords;

    // the row order has to match the node and shape-function ordering
    local_coords << Precision(0), Precision(0),  // node 1
                    Precision(1), Precision(0),  // node 2
                    Precision(0), Precision(1);  // node 3

    return local_coords;
}

/**
 * @brief Computes the closest point on the triangular boundary.
 *
 * The three triangle edges are represented by linear boundary elements. The
 * global point is projected onto every edge, after which the closest projection
 * is transformed back into the natural coordinates of the triangle.
 *
 * @param global Global point to project onto the element boundary.
 * @param node_coords Global coordinates of the three triangle nodes.
 *
 * @return Natural triangle coordinates of the closest boundary point.
 */
Vec2 Surface3::closest_point_on_boundary(const Vec3&               global,
                                         const StaticMatrix<3, 3>& node_coords) const {
    // represent the three triangle edges by line elements following the
    // surface node ordering
    Line2B edge_01({0, 1});
    Line2B edge_12({1, 2});
    Line2B edge_20({2, 0});

    // the line-element interface operates on Field storage, so transfer the
    // supplied fixed-size coordinate matrix into a temporary nodal field
    Field node_field("SURFACE3_BOUNDARY", FieldDomain::NODE, 3, 3);

    for (Index local_id = 0; local_id < 3; ++local_id) {
        for (Dim component = 0; component < 3; ++component) {
            node_field(local_id, component) = node_coords(local_id, component);
        }
    }

    // project the global point independently onto every triangle edge
    const Precision edge_01_local = edge_01.global_to_local(global, node_field);
    const Precision edge_12_local = edge_12.global_to_local(global, node_field);
    const Precision edge_20_local = edge_20.global_to_local(global, node_field);

    // map the edge-local projections back into physical coordinates so their
    // distances to the requested global point can be compared
    const Vec3 point_01 = edge_01.local_to_global(edge_01_local, node_field);
    const Vec3 point_12 = edge_12.local_to_global(edge_12_local, node_field);
    const Vec3 point_20 = edge_20.local_to_global(edge_20_local, node_field);

    // squared distances are sufficient for comparison and avoid unnecessary
    // square-root evaluations
    const Precision distance_01 = (point_01 - global).squaredNorm();
    const Precision distance_12 = (point_12 - global).squaredNorm();
    const Precision distance_20 = (point_20 - global).squaredNorm();

    // edge 0-1 follows (r,s) = (p,0)
    if (distance_01 <= distance_12 && distance_01 <= distance_20) {
        return {edge_01_local, Precision(0)};
    }

    // edge 2-0 follows (r,s) = (0,1-p) because its line orientation starts at
    // node 2 and ends at node 0
    if (distance_20 <= distance_01 && distance_20 <= distance_12) {
        return {Precision(0), Precision(1) - edge_20_local};
    }

    // edge 1-2 follows (r,s) = (1-p,p)
    return {
        Precision(1) - edge_12_local,
        edge_12_local
    };
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
bool Surface3::in_bounds(const Vec2& local) const {
    const Precision r = local(0);
    const Precision s = local(1);

    return r >= Precision(0) &&
           s >= Precision(0) &&
           r + s <= Precision(1);
}

/**
 * @brief Returns the quadrature rule used for the linear triangle.
 *
 * A linear rule is sufficient for the interpolation order of the element.
 * The static object is constructed only once and reused for all evaluations.
 *
 * @return Quadrature rule on the isoparametric triangular domain.
 */
const fem::math::quadrature::Quadrature& Surface3::integration_scheme() const {
    static const fem::math::quadrature::Quadrature scheme{
        fem::math::quadrature::DOMAIN_ISO_TRI,
        fem::math::quadrature::ORDER_LINEAR
    };

    return scheme;
}

} // namespace fem::model