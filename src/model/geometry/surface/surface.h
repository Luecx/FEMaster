/**
 * @file surface.h
 * @brief Defines the common base template for finite-element surfaces.
 *
 * The surface abstraction stores the connectivity of one triangular or
 * quadrilateral finite element and provides the topology-independent
 * operations needed by geometry, projection and surface integration.
 * Concrete surface types supply the shape functions, their derivatives, the
 * natural nodal coordinates, the boundary projection and the quadrature rule.
 *
 * The implementation covers interpolation, mappings between natural and
 * global coordinates, differential geometry, polygon clipping and integration
 * of scalar, vector and tensor fields over the physical surface.
 *
 * @see SurfaceInterface
 * @see SurfacePolygon
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "surface_interface.h"

#include <array>
#include <functional>

#include <Eigen/Geometry>

namespace fem::model {

/**
 * @brief Common base class for triangular and quadrilateral surface elements.
 *
 * `N = 3` and `N = 6` represent linear and quadratic triangles. `N = 4` and
 * `N = 8` represent bilinear and serendipity quadrilaterals. The natural
 * coordinate domain is therefore topology-dependent, while all physical
 * geometry is obtained from the isoparametric mapping defined by the derived
 * shape functions.
 *
 * The class owns the global node identifiers and uses them to gather nodal
 * coordinates from a model field. Derived classes must implement the
 * element-specific interpolation data and the closest-point search on the
 * element boundary.
 *
 * @tparam N Number of surface nodes.
 */
template<Index N>
struct Surface : public SurfaceInterface {
    static_assert(N == 3 || N == 4 || N == 6 || N == 8,
        "Surface supports only 3-, 4-, 6- and 8-node elements");

    static constexpr Index num_edges          = N > 4 ? N / 2 : N;
    static constexpr Index num_nodes          = N;
    static constexpr Index num_nodes_per_edge = N > 4 ? 3 : 2;

    // Global node identifiers of the surface
    std::array<ID, N> nodeIds{};

    // Construction
    explicit Surface(const std::array<ID, N>& node_ids);

    ~Surface() override = default;

    // Shape functions
    virtual StaticMatrix<N, 1> shape_function         (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> shape_derivative       (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 3> shape_second_derivative(Precision r, Precision s) const = 0;

    // Evaluate the shape functions through the dynamic interface
    DynamicVector shape_function(const Vec2& local) const override;

    // Natural coordinates of the surface nodes
    virtual StaticMatrix<N, 2> node_coords_local() const = 0;

    // Gather global coordinates of the attached nodes
    StaticMatrix<N, 3> node_coords_global(const Field& node_coords) const;

    // Return contiguous mutable connectivity storage
    ID* nodes() override;

    // Boundary projection delegated to the concrete surface implementation
    virtual Vec2 closest_point_on_boundary(
        const Vec3&               global,
        const StaticMatrix<N, 3>& node_coords) const = 0;

    // Numerical integration scheme
    virtual const math::quadrature::Quadrature& integration_scheme() const = 0;

    // Interpolation and differential geometry
    template<int M>
    StaticVector<M> interpolate(const StaticMatrix<N, M>& nodal_values, Precision r, Precision s) const;
    StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords , Precision r, Precision s) const;

    // Return the polygon of the element domain in natural coordinates
    Polygon local_domain_polygon() const override;

    // Coordinate transformations
    Vec3 local_to_global(const Vec2& local, const Field& node_coords) const override;
    Vec2 global_to_local(const Vec3& global, const Field& node_coords, bool clip = false) const override;

    // Surface geometry
    Vec3      normal(const Field& node_coords, const Vec2& local) const override;
    Precision area  (const Field& node_coords) const override;

    // Surface integration
    Precision integrate_scalar_field(
        const Field&       node_coords,
        const ScalarField& field) const override;
    Vec3 integrate_vector_field(
        const Field&    node_coords,
        const VecField& field) const override;
    void integrate_vector_field(
        const Field&    node_coords,
              Field&    target,
        const VecField& field) const override;
    Mat3 integrate_tensor_field(
        const Field&    node_coords,
        const TenField& field) const override;
    void integrate_triangular(
        const Field&                                                    node_coords,
        const Polygon&                                                  polygon,
        const math::quadrature::Quadrature&                             scheme,
        const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand
    ) const override;
};

} // namespace fem::model

#include "surface_geometry.inl"
#include "surface_projection.inl"
#include "surface_integrate.inl"
