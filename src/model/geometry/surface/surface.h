/**
 * @file surface.h
 * @brief Defines the common base template for finite-element surfaces.
 *
 * The surface abstraction stores the connectivity of one triangular or
 * quadrilateral finite-element surface and provides the topology-independent
 * operations needed by geometry, projection and numerical integration. Concrete
 * surface types supply the shape functions, their derivatives, natural nodal
 * coordinates, boundary projection and the ordinary surface quadrature rule.
 *
 * The common implementation covers interpolation, mappings between natural and
 * global coordinates, differential geometry, polygon clipping and integration
 * of scalar, vector and tensor fields over the physical surface. In addition to
 * ordinary field integration, scalar boundary formulations may assemble
 * consistent nodal vectors or scalar-weighted shape-function product matrices.
 * The latter uses a dedicated quadrature rule so higher-order `N_i N_j` products
 * do not alter the established integration behavior of mechanical surface loads.
 *
 * Surface geometry is used by structural pressure and traction loads, thermal
 * heat-flux and convection conditions, contact, coupling and other algorithms
 * operating on compiled element boundaries.
 *
 * @see SurfaceInterface
 * @see Polygon
 * @see surface_integrate.inl
 * @see surface_shape_matrix.inl
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
 * element-specific interpolation data, ordinary quadrature rule and closest-
 * point search on the element boundary. This base class implements the common
 * coordinate mappings, physical Jacobians, normals, areas and surface-field
 * integration routines.
 *
 * Scalar integration is available in three related forms: direct integration of
 * a scalar field, consistent distribution into a one-component nodal field and
 * integration of `N N^T` weighted by a scalar field. These operations share the
 * same physical surface geometry but deliberately need not share quadrature
 * order.
 *
 * @tparam N Number of surface nodes.
 */
template<Index N>
struct Surface : public SurfaceInterface {
    static_assert(N == 3 || N == 4 || N == 6 || N == 8,
        "Surface supports only 3-, 4-, 6- and 8-node elements");

    // Compile-time topology information used by boundary extraction and generic
    // surface algorithms. Quadratic surfaces carry one midside node per edge.
    static constexpr Index num_edges          = N > 4 ? N / 2 : N;
    static constexpr Index num_nodes          = N;
    static constexpr Index num_nodes_per_edge = N > 4 ? 3 : 2;

    // Global node identifiers of the compiled surface. Their ordering must match
    // the concrete shape-function vector and natural nodal coordinates.
    std::array<ID, N> nodeIds{};

    // Construction from compiled connectivity
    explicit Surface(const std::array<ID, N>& node_ids);
    ~Surface() override = default;

    // Shape functions and their first and second derivatives in natural
    // coordinates. Concrete surface types provide the topology-specific
    // interpolation while the base class consumes it for all common operations.
    virtual StaticMatrix<N, 1> shape_function         (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> shape_derivative       (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 3> shape_second_derivative(Precision r, Precision s) const = 0;

    // Dynamic-interface views of the interpolation data. These retain the same
    // node ordering as the fixed-size concrete shape-function implementation.
    DynamicVector shape_function(const Vec2& local) const override;
    DynamicMatrix node_coords_natural() const override;

    // Natural and global nodal coordinates. The natural ordering must match
    // `nodeIds`; global coordinates are gathered from the supplied NODE field.
    virtual StaticMatrix<N, 2> node_coords_local() const = 0;
    StaticMatrix<N, 3> node_coords_global(const Field& node_coords) const;

    // Contiguous connectivity access used by topology-independent model code.
    ID* nodes() override;
    const ID* nodes() const override;

    // Project a physical point onto the element boundary. The concrete topology
    // owns the edge parameterization while the base class owns interior mapping.
    virtual Vec2 closest_point_on_boundary(
        const Vec3& global,
        const StaticMatrix<N, 3>& node_coords) const = 0;

    // Ordinary quadrature used for standard surface-field integration. Dedicated
    // shape-product integration selects its own higher-order rule separately.
    virtual const math::quadrature::Quadrature& integration_scheme() const = 0;

    // Isoparametric interpolation and physical differential geometry. The
    // Jacobian columns are the physical tangent vectors associated with the two
    // natural coordinates.
    template<int M>
    StaticVector<M> interpolate(const StaticMatrix<N, M>& nodal_values, Precision r, Precision s) const;
    StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords, Precision r, Precision s) const;

    // Polygon representing the valid natural element domain. It is used by
    // clipping and polygon-restricted integration routines.
    Polygon local_domain_polygon() const override;

    // Coordinate transformations between natural surface coordinates and the
    // global physical surface. Global-to-local inversion may optionally clip the
    // recovered point to the valid element domain.
    Vec3 local_to_global(const Vec2& local, const Field& node_coords) const override;
    Vec2 global_to_local(const Vec3& global, const Field& node_coords, bool clip = false) const override;

    // Surface geometry in the physical configuration represented by
    // `node_coords`. The normal follows the connectivity orientation.
    Vec3 normal(const Field& node_coords, const Vec2& local) const override;
    Precision area(const Field& node_coords) const override;

    // Surface integration of physical fields. Scalar and vector overloads with a
    // target field perform consistent nodal distribution using the surface shape
    // functions; direct-return overloads integrate only the supplied field value.
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

    // Integrate a scalar-weighted `N N^T` matrix. This operation is used by
    // Robin boundary operators and employs a dedicated quadrature order suitable
    // for products of interpolation functions.
    DynamicMatrix integrate_scalar_shape_matrix(
        const Field& node_coords,
        const ScalarField& field) const override;

    // Integrate over the overlap of a supplied natural-coordinate polygon with
    // the valid element domain. The callback receives natural coordinates,
    // physical coordinates and the complete physical integration weight.
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
#include "surface_shape_matrix.inl"
