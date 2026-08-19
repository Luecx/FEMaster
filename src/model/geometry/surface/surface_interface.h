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

struct SurfaceInterface {
    using Ptr = std::shared_ptr<SurfaceInterface>;

    using ScalarField = ::fem::ScalarField;
    using VecField    = ::fem::VecField;
    using TenField    = ::fem::TenField;
    using Polygon     = SurfacePolygon<4>;

    const Index n_edges;
    const Index n_nodes;
    const Index n_nodes_per_edge;

    SurfaceInterface(Index edge_count = 0,
                     Index node_count = 0,
                     Index edge_node_count = 0)
        : n_edges(edge_count),
          n_nodes(node_count),
          n_nodes_per_edge(edge_node_count) {}

    virtual ~SurfaceInterface() = default;

    virtual ID* nodes() = 0;
    virtual const ID* nodes() const = 0;
    virtual Ptr copy() const = 0;

    ID* begin() { return nodes(); }
    ID* end() { return nodes() + n_nodes; }
    const ID* begin() const { return nodes(); }
    const ID* end() const { return nodes() + n_nodes; }

    virtual Vec3 local_to_global(const Vec2& local, const Field& node_coords) const = 0;
    virtual Vec2 global_to_local(const Vec3& global, const Field& node_coords, bool clip = false) const = 0;

    virtual Vec3 normal(const Field& node_coords, const Vec2& local) const = 0;
    virtual bool in_bounds(const Vec2& local) const = 0;
    virtual Precision area(const Field& node_coords) const = 0;

    virtual Polygon local_domain_polygon() const = 0;
    virtual DynamicVector shape_function(const Vec2& local) const = 0;
    virtual DynamicMatrix node_coords_natural() const = 0;

    virtual Precision integrate_scalar_field(const Field& node_coords, const ScalarField& field) const = 0;
    virtual Vec3 integrate_vector_field(const Field& node_coords, const VecField& field) const = 0;
    virtual void integrate_vector_field(const Field& node_coords, Field& target, const VecField& field) const = 0;
    virtual Mat3 integrate_tensor_field(const Field& node_coords, const TenField& field) const = 0;
    virtual void integrate_triangular(
        const Field& node_coords,
        const Polygon& polygon,
        const math::quadrature::Quadrature& scheme,
        const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand) const = 0;
};

} // namespace fem::model
