/**
 * @file surface.h
 * @brief Defines the common base template for finite-element surfaces.
 */

#pragma once

#include "surface_interface.h"

#include <array>
#include <functional>

#include <Eigen/Geometry>

namespace fem::model {

template<Index N>
struct Surface : public SurfaceInterface {
    static_assert(N == 3 || N == 4 || N == 6 || N == 8,
        "Surface supports only 3-, 4-, 6- and 8-node elements");

    static constexpr Index num_edges          = N > 4 ? N / 2 : N;
    static constexpr Index num_nodes          = N;
    static constexpr Index num_nodes_per_edge = N > 4 ? 3 : 2;

    std::array<ID, N> nodeIds{};

    explicit Surface(const std::array<ID, N>& node_ids);
    ~Surface() override = default;

    virtual StaticMatrix<N, 1> shape_function         (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> shape_derivative       (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 3> shape_second_derivative(Precision r, Precision s) const = 0;

    DynamicVector shape_function(const Vec2& local) const override;
    DynamicMatrix node_coords_natural() const override;

    virtual StaticMatrix<N, 2> node_coords_local() const = 0;
    StaticMatrix<N, 3> node_coords_global(const Field& node_coords) const;

    ID* nodes() override;
    const ID* nodes() const override;

    virtual Vec2 closest_point_on_boundary(
        const Vec3& global,
        const StaticMatrix<N, 3>& node_coords) const = 0;

    virtual const math::quadrature::Quadrature& integration_scheme() const = 0;

    template<int M>
    StaticVector<M> interpolate(const StaticMatrix<N, M>& nodal_values, Precision r, Precision s) const;
    StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords, Precision r, Precision s) const;

    Polygon local_domain_polygon() const override;

    Vec3 local_to_global(const Vec2& local, const Field& node_coords) const override;
    Vec2 global_to_local(const Vec3& global, const Field& node_coords, bool clip = false) const override;

    Vec3 normal(const Field& node_coords, const Vec2& local) const override;
    Precision area(const Field& node_coords) const override;

    Precision integrate_scalar_field(
        const Field& node_coords,
        const ScalarField& field) const override;
    void integrate_scalar_field(
        const Field& node_coords,
        Field& target,
        const ScalarField& field) const override;
    Vec3 integrate_vector_field(
        const Field& node_coords,
        const VecField& field) const override;
    void integrate_vector_field(
        const Field& node_coords,
        Field& target,
        const VecField& field) const override;
    Mat3 integrate_tensor_field(
        const Field& node_coords,
        const TenField& field) const override;
    DynamicMatrix integrate_scalar_shape_matrix(
        const Field& node_coords,
        const ScalarField& field) const override;
    void integrate_triangular(
        const Field& node_coords,
        const Polygon& polygon,
        const math::quadrature::Quadrature& scheme,
        const std::function<void(const Vec2&, const Vec3&, Precision)>& integrand
    ) const override;
};

} // namespace fem::model

#include "surface_geometry.inl"
#include "surface_projection.inl"
#include "surface_integrate.inl"
