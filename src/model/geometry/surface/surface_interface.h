
/**
 * @file surface_interface.h
 * @brief Declares the common interface for surface geometry, integration and loading.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "../../../core/core.h"
#include "../../../data/field.h"
#include "../../../math/quadrature.h"

#include <array>
#include <memory>
#include <vector>

namespace fem::model {

/**
 * @brief Common interface for finite-element surface geometry and loading operations.
 */
struct SurfaceInterface {
    using Ptr = std::shared_ptr<SurfaceInterface>;

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
    virtual Vec3              normal              (const Field& node_coords, const Vec2& local) const = 0;
    virtual bool              in_bounds           (const Vec2& local) const                            = 0;
    virtual Precision         area                (const Field& node_coords) const                     = 0;
    virtual Precision         jacobian_measure    (const Field& node_coords, const Vec2& local) const = 0;
    virtual std::vector<Vec2> local_domain_polygon() const                                             = 0;

    // Shape functions
    virtual DynamicVector shape_function         (const Vec2& local) const             = 0;
    virtual DynamicVector shape_function_integral(const Field& node_coords) const      = 0;

    // Connectivity
    virtual ID* nodes() = 0;

    // Surface loads
    virtual void apply_dload(const Field& node_coords,
                             Field&       node_loads,
                             Vec3         load) = 0;

    virtual void apply_pload(const Field& node_coords,
                             Field&       node_loads,
                             Precision    load) = 0;

    /**
     * @brief Integrates the portion of a local triangle inside the surface domain.
     *
     * The triangle is clipped against the local element domain, triangulated and
     * integrated using the supplied quadrature scheme.
     *
     * @param triangle   Triangle vertices in local coordinates.
     * @param node_coords Global nodal coordinates.
     * @param scheme     Quadrature scheme for the resulting triangles.
     * @param integrand  Callable receiving local position, global position and weight.
     */
    template<typename F>
    void integrate_local_triangle(const std::array<Vec2, 3>&    triangle,
                                  const Field&                   node_coords,
                                  const quadrature::Quadrature& scheme,
                                  const F&                       integrand) const;

    ID* begin() { return nodes(); }
    ID* end()   { return nodes() + n_nodes; }
};

} // namespace fem::model

#include "surface_interface.inl"

