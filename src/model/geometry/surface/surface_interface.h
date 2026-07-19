/**
 * @file surface_interface.h
 * @brief Declares the common interface for surface geometry, integration and loading.
 *
 * The interface defines the operations required by triangular and
 * quadrilateral finite-element surfaces independently of their interpolation
 * order. It covers mappings between natural and global coordinates, domain
 * checks, surface geometry, field integration and connectivity access.
 * Concrete surface types provide the element-specific shape functions,
 * boundary behavior and integration implementation.
 *
 * @see Surface
 * @see SurfacePolygon
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

#include <array>
#include <memory>
#include <vector>

namespace fem::model {

/**
 * @brief Common interface for finite-element surface geometry and loading operations.
 *
 * A surface implementation exposes its topology through the edge and node
 * counts and uses a shared set of natural-coordinate and global-coordinate
 * operations. Surface fields are evaluated in global coordinates, while
 * element quadrature is defined in the natural domain and mapped to physical
 * surface measure by the concrete implementation.
 *
 * Derived classes must provide the coordinate transformations, geometric
 * queries, field integration routines, natural-domain polygon and contiguous
 * node connectivity expected by the model and solver infrastructure.
 */
struct SurfaceInterface {
    using Ptr = std::shared_ptr<SurfaceInterface>;

    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;

    // Natural-domain polygons use only the element corner coordinates. Four
    // vertices are therefore sufficient for both triangular and quadrilateral
    // surfaces, including polygons supplied for restricted integration.
    using Polygon = SurfacePolygon<4>;

    // Topological counts used by connectivity, boundary projection and generic
    // surface algorithms.
    const Index n_edges;
    const Index n_nodes;
    const Index n_nodes_per_edge;

    SurfaceInterface(Index edge_count      = 0,
                     Index node_count      = 0,
                     Index edge_node_count = 0)
        : n_edges         (edge_count),
          n_nodes         (node_count),
          n_nodes_per_edge(edge_node_count) {}

    virtual ~SurfaceInterface() = default;

    // Contiguous global node connectivity used by generic model traversal and
    // range-based access through begin() and end().
    virtual ID* nodes() = 0;

    ID* begin() { return nodes(); }
    ID* end()   { return nodes() + n_nodes; }

    // Coordinate transformations between the natural element domain and the
    // physical global coordinates represented by the nodal field.
    virtual Vec3 local_to_global(const Vec2& local,
                                 const Field& node_coords) const = 0;

    virtual Vec2 global_to_local(const Vec3& global,
                                 const Field& node_coords,
                                 bool         clip = false) const = 0;

    // Surface geometry queries. `normal` is evaluated at a natural coordinate,
    // while `area` integrates the physical measure over the complete element.
    virtual Vec3       normal   (const Field& node_coords, const Vec2& local) const  = 0;
    virtual bool       in_bounds(const Vec2& local) const                            = 0;
    virtual Precision  area     (const Field& node_coords) const                     = 0;

    // Return the valid natural-coordinate domain as a counter-clockwise
    // polygon. Triangles use the reference triangle; quadrilaterals use the
    // reference square. The polygon is also used for clipping restricted
    // integration regions.
    virtual Polygon local_domain_polygon() const = 0;

    // Evaluate the element shape functions at one natural-coordinate point.
    // The dynamic return type allows generic surface algorithms to consume
    // different interpolation orders through the common interface.
    virtual DynamicVector shape_function         (const Vec2& local) const             = 0;

    // Integrate scalar, vector and tensor fields over the physical surface.
    // Field callbacks receive interpolated global positions, and the concrete
    // implementation applies the physical surface Jacobian to the quadrature
    // measure.
    virtual Precision integrate_scalar_field(
        const Field&       node_coords,
        const ScalarField& field) const = 0;
    virtual Vec3 integrate_vector_field(
        const Field&    node_coords,
        const VecField& field) const = 0;
    virtual void integrate_vector_field(
        const Field&    node_coords,
        Field&          target,
        const VecField& field) const = 0;
    virtual Mat3 integrate_tensor_field(
        const Field& node_coords,
        const TenField&    field) const = 0;
    virtual void integrate_triangular(
        const Field& node_coords,
        const Polygon& polygon,
        const math::quadrature::Quadrature& scheme,
        const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand) const = 0;

};

} // namespace fem::model
