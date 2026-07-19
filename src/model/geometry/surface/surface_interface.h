/**
 * @file surface_interface.h
 * @brief Declares the common interface for surface geometry, integration and loading.
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
 */
struct SurfaceInterface {
    using Ptr = std::shared_ptr<SurfaceInterface>;

    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;

    // define every polygon to have max 4 points, this is sufficient as for quads we only use corner points
    // as its in local coordinates. Whatever we provide as the polygon for the integration must also be 4 nodes max.
    using Polygon = SurfacePolygon<4>;

    // Number of edges
    const Index n_edges;

    // Number of nodes
    const Index n_nodes;

    // Number of nodes per edge
    const Index n_nodes_per_edge;

    SurfaceInterface(Index edge_count      = 0,
                     Index node_count      = 0,
                     Index edge_node_count = 0)
        : n_edges         (edge_count),
          n_nodes         (node_count),
          n_nodes_per_edge(edge_node_count) {}

    virtual ~SurfaceInterface() = default;

    // Coordinate transformations
    virtual Vec3 local_to_global(const Vec2& local,
                                 const Field& node_coords) const = 0;

    virtual Vec2 global_to_local(const Vec3& global,
                                 const Field& node_coords,
                                 bool         clip = false) const = 0;

    // Surface geometry
    virtual Vec3       normal              (const Field& node_coords, const Vec2& local) const  = 0;
    virtual bool       in_bounds           (const Vec2& local) const                            = 0;
    virtual Precision  area                (const Field& node_coords) const                     = 0;

    // polygon of the local domain which is either the [-1,1]x[-1,1]
    // or the triangular region in natural coordinates
    virtual Polygon local_domain_polygon() const = 0;

    // Shape functions
    virtual DynamicVector shape_function         (const Vec2& local) const             = 0;

    // integrating over the surface
    virtual Precision integrate_scalar_field(const Field& node_coords, const ScalarField& field) const = 0;
    virtual Vec3      integrate_vector_field(const Field& node_coords, const VecField&    field) const = 0;
    virtual void      integrate_vector_field(const Field& node_coords,       Field&       target, const VecField& field) const = 0;
    virtual Mat3      integrate_tensor_field(const Field& node_coords, const TenField&    field) const = 0;
    virtual void      integrate_triangular  (const Field& node_coords,
                                             const Polygon& polygon,
                                             const math::quadrature::Quadrature& scheme,
                                             const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand) const = 0;

    // Connectivity
    virtual ID* nodes() = 0;

    ID* begin() { return nodes(); }
    ID* end()   { return nodes() + n_nodes; }
};

} // namespace fem::model
