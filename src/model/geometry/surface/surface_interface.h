/**
 * @file surface_interface.h
 * @brief Declares the topology-independent interface for finite-element surfaces.
 *
 * `SurfaceInterface` allows model algorithms to operate on triangular and
 * quadrilateral surfaces of different interpolation orders without knowing their
 * compile-time node count. It exposes connectivity, coordinate transformations,
 * physical geometry, shape interpolation and the integration operations needed
 * by surface loads, thermal boundary conditions, contact and coupling.
 *
 * Concrete fixed-size `Surface<N>` implementations provide the numerical
 * algorithms. This interface only standardizes their dynamic polymorphic view,
 * including scalar/vector/tensor integration, consistent nodal distribution and
 * scalar-weighted `N N^T` matrix integration for Robin boundary operators.
 *
 * @see Surface
 * @see SurfacePolygon
 * @see surface_integrate.inl
 * @see surface_shape_matrix.inl
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "../../../core/core.h"
#include "../../../core/types_eig.h"
#include "../../../data/field.h"
#include "../../../math/quadrature.h"
#include "surface_polygon.h"

#include <functional>
#include <memory>

namespace fem::model {

/**
 * @brief Dynamic interface shared by all finite-element surface topologies.
 *
 * A surface exposes the global identifiers of its compiled nodes together with
 * the geometry and integration operations defined by its interpolation. The
 * three immutable topology counts describe the number of edges, total nodes and
 * nodes per edge, allowing callers to traverse connectivity without downcasting
 * to a particular `Surface<N>` specialization.
 *
 * Coordinate mappings use natural two-dimensional surface coordinates and
 * physical three-dimensional global coordinates. Geometry methods evaluate the
 * physical surface represented by the nodal coordinate field supplied by the
 * caller; the same interface can therefore be used with reference or current
 * positions when an algorithm explicitly selects the desired configuration.
 *
 * Integration callbacks are evaluated in global coordinates. Overloads with a
 * target `Field` perform consistent nodal distribution through the surface shape
 * functions, whereas direct-return overloads integrate the supplied field itself.
 * `integrate_scalar_shape_matrix()` constructs a dense local scalar operator and
 * is intended for terms such as Robin convection matrices.
 */
struct SurfaceInterface {
    // Shared ownership type used by compiled model surface storage and boundary-
    // condition regions containing heterogeneous surface topologies.
    using Ptr = std::shared_ptr<SurfaceInterface>;

    // Physical callback types used by generic surface integration. Functions are
    // evaluated at global quadrature-point positions.
    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;

    // Polygon type used for natural-domain clipping and polygon-restricted
    // integration. Four vertices are sufficient for the supported complete
    // triangular and quadrilateral natural domains.
    using Polygon = SurfacePolygon<4>;

    // Immutable topology counts supplied by the concrete fixed-size surface at
    // construction. Quadratic surfaces carry three nodes per edge.
    const Index n_edges;
    const Index n_nodes;
    const Index n_nodes_per_edge;

    // Construction from topology metadata. The interface itself owns no node
    // identifiers; concrete surfaces provide contiguous connectivity storage.
    SurfaceInterface(Index edge_count = 0,
                     Index node_count = 0,
                     Index edge_node_count = 0)
        : n_edges(edge_count),
          n_nodes(node_count),
          n_nodes_per_edge(edge_node_count) {}

    virtual ~SurfaceInterface() = default;

    // Connectivity access and polymorphic copying. The returned node storage
    // contains exactly `n_nodes` global compiled identifiers in shape-function
    // ordering.
    virtual ID* nodes() = 0;
    virtual const ID* nodes() const = 0;
    virtual Ptr copy() const = 0;

    // Iteration over contiguous surface connectivity. These helpers intentionally
    // expose the same ordering as `nodes()` and the surface shape-function vector.
    ID* begin() { return nodes(); }
    ID* end() { return nodes() + n_nodes; }
    const ID* begin() const { return nodes(); }
    const ID* end() const { return nodes() + n_nodes; }

    // Coordinate transformations between natural element coordinates and the
    // physical surface defined by `node_coords`. Global-to-local inversion may
    // optionally clip the result to the valid natural element domain.
    virtual Vec3 local_to_global(const Vec2& local, const Field& node_coords) const = 0;
    virtual Vec2 global_to_local(const Vec3& global, const Field& node_coords, bool clip = false) const = 0;

    // Physical surface geometry. The normal follows connectivity orientation;
    // `in_bounds` evaluates the natural domain and `area` integrates the physical
    // surface measure represented by the supplied nodal coordinates.
    virtual Vec3 normal(const Field& node_coords, const Vec2& local) const = 0;
    virtual bool in_bounds(const Vec2& local) const = 0;
    virtual Precision area(const Field& node_coords) const = 0;

    // Natural-domain and interpolation data in dynamic form. Natural nodal
    // coordinates and shape functions use exactly the same node ordering as the
    // global connectivity returned by `nodes()`.
    virtual Polygon local_domain_polygon() const = 0;
    virtual DynamicVector shape_function(const Vec2& local) const = 0;
    virtual DynamicMatrix node_coords_natural() const = 0;

    // Physical surface-field integration. Direct-return overloads integrate the
    // supplied scalar/vector/tensor quantity; target-field overloads multiply by
    // surface shape functions and scatter consistent nodal contributions.
    virtual Precision integrate_scalar_field(const Field& node_coords, const ScalarField& field) const = 0;
    virtual void integrate_scalar_field(const Field& node_coords, Field& target, const ScalarField& field) const = 0;
    virtual Vec3 integrate_vector_field(const Field& node_coords, const VecField& field) const = 0;
    virtual void integrate_vector_field(const Field& node_coords, Field& target, const VecField& field) const = 0;
    virtual Mat3 integrate_tensor_field(const Field& node_coords, const TenField& field) const = 0;

    // Scalar-weighted shape-product integration for local boundary operators of
    // the form integral a(x) N N^T dGamma. Concrete surfaces choose a quadrature
    // rule appropriate for the higher polynomial order of the shape products.
    virtual DynamicMatrix integrate_scalar_shape_matrix(
        const Field& node_coords,
        const ScalarField& field) const = 0;

    // Polygon-restricted integration in natural coordinates. Implementations clip
    // the supplied polygon against the valid element domain and pass natural
    // coordinates, physical coordinates and the complete area weight to the
    // integrand callback.
    virtual void integrate_triangular(
        const Field& node_coords,
        const Polygon& polygon,
        const math::quadrature::Quadrature& scheme,
        const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand) const = 0;
};

} // namespace fem::model
