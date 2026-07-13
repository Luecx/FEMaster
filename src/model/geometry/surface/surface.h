
/**
 * @file surface.h
 * @brief Defines the common base template for finite-element surfaces.
 *
 * @author Finn Eggers
 * @date 27.09.2024
 */

#pragma once

#include "surface_interface.h"

#include <Eigen/Cholesky>

#include <array>
#include <functional>
#include <limits>
#include <vector>

namespace fem::model {

/**
 * @brief Common base class for triangular and quadrilateral surface elements.
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

    // Global node identifiers
    std::array<ID, N> nodeIds{};

    explicit Surface(const std::array<ID, N>& node_ids)
        : SurfaceInterface(num_edges, num_nodes, num_nodes_per_edge),
          nodeIds         (node_ids) {}

    ~Surface() override = default;

    // Shape functions and local geometry
    virtual StaticMatrix<N, 1> shape_function         (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> shape_derivative       (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 3> shape_second_derivative(Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> node_coords_local      () const                          = 0;

    // Boundary projection
    virtual Vec2 closest_point_on_boundary(const Vec3&               global,
                                           const StaticMatrix<N, 3>& node_coords) const = 0;

    // Numerical integration
    virtual const quadrature::Quadrature& integration_scheme() const = 0;

    /**
     * @brief Interpolates nodal values at the local coordinates `(r, s)`.
     */
    template<Index M>
    StaticVector<M> interpolate(const StaticMatrix<N, M>& nodal_values,
                                Precision                 r,
                                Precision                 s) const {
        return (shape_function(r, s).transpose() * nodal_values).transpose();
    }

    /**
     * @brief Evaluates the shape functions using dynamic storage.
     */
    DynamicVector shape_function(const Vec2& local) const override {
        return DynamicVector(shape_function(local(0), local(1)));
    }

    /**
     * @brief Computes the local-to-global surface Jacobian.
     */
    StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords,
                                Precision                 r,
                                Precision                 s) const {
        return node_coords.transpose() * shape_derivative(r, s);
    }

    /**
     * @brief Computes the local surface area scaling.
     */
    Precision jacobian_measure(const Field& node_coords,
                               const Vec2&  local) const override {
        const auto coordinates = node_coords_global(node_coords);
        const auto jac         = jacobian(coordinates, local(0), local(1));

        return jac.col(0).cross(jac.col(1)).norm();
    }

    /**
     * @brief Returns the polygon of the local element domain.
     */
    std::vector<Vec2> local_domain_polygon() const override {
        if constexpr (N == 3 || N == 6) {
            return {
                Vec2(Precision(0), Precision(0)),
                Vec2(Precision(1), Precision(0)),
                Vec2(Precision(0), Precision(1))
            };
        } else {
            return {
                Vec2(Precision(-1), Precision(-1)),
                Vec2(Precision( 1), Precision(-1)),
                Vec2(Precision( 1), Precision( 1)),
                Vec2(Precision(-1), Precision( 1))
            };
        }
    }

    /**
     * @brief Maps local coordinates onto the global surface.
     */
    Vec3 local_to_global(const Vec2& local,
                         const Field& node_coords) const override {
        const auto coordinates = node_coords_global(node_coords);
        return coordinates.transpose() * shape_function(local(0), local(1));
    }

    /**
     * @brief Projects a global point onto the surface.
     *
     * Newton iterations are started from several local coordinates. When clipping
     * is enabled, only interior projections and the closest boundary projection
     * are considered.
     */
    Vec2 global_to_local(const Vec3& global,
                         const Field& node_coords,
                         bool         clip = false) const override {
        constexpr Index     max_iterations = 32;
        constexpr Precision tolerance      = Precision(1e-12);

        // Select initial guesses for the element topology
        std::vector<Vec2> initial_guesses{};

        if constexpr (N == 3) {
            initial_guesses = {
                {Precision(0.25), Precision(0.25)}
            };
        } else if constexpr (N == 4) {
            initial_guesses = {
                {Precision(-1), Precision(-1)},
                {Precision( 1), Precision(-1)},
                {Precision( 1), Precision( 1)},
                {Precision(-1), Precision( 1)}
            };
        } else if constexpr (N == 6) {
            initial_guesses = {
                {Precision(0), Precision(0)},
                {Precision(1), Precision(0)},
                {Precision(0), Precision(1)}
            };
        } else {
            initial_guesses = {
                {Precision(-1), Precision(-1)},
                {Precision( 1), Precision(-1)},
                {Precision( 1), Precision( 1)},
                {Precision(-1), Precision( 1)},
                {Precision( 0), Precision(-1)},
                {Precision( 1), Precision( 0)},
                {Precision( 0), Precision( 1)},
                {Precision(-1), Precision( 0)}
            };
        }

        const auto coordinates = node_coords_global(node_coords);

        auto squared_distance = [&](const Vec2& local) {
            const Vec3 position =
                coordinates.transpose() * shape_function(local(0), local(1));

            return (position - global).squaredNorm();
        };

        Vec2      best_local            = initial_guesses.front();
        Precision best_distance_squared = std::numeric_limits<Precision>::max();

        // Project each initial guess with Newton iterations
        for (const Vec2& initial_guess : initial_guesses) {
            Precision r = initial_guess(0);
            Precision s = initial_guess(1);

            for (Index iteration = 0; iteration < max_iterations; ++iteration) {
                const auto shape             = shape_function(r, s);
                const auto derivative        = shape_derivative(r, s);
                const auto second_derivative = shape_second_derivative(r, s);

                const Vec3               position = coordinates.transpose() * shape;
                const StaticMatrix<3, 2> first    = coordinates.transpose() * derivative;
                const StaticMatrix<3, 3> second   = coordinates.transpose() * second_derivative;

                const Vec3 dx_dr    = first.col(0);
                const Vec3 dx_ds    = first.col(1);
                const Vec3 d2x_dr2  = second.col(0);
                const Vec3 d2x_ds2  = second.col(1);
                const Vec3 d2x_drds = second.col(2);

                const Vec3 difference = position - global;

                const Precision gradient_r = difference.dot(dx_dr);
                const Precision gradient_s = difference.dot(dx_ds);

                const Precision hessian_rr = dx_dr.dot(dx_dr) + difference.dot(d2x_dr2);
                const Precision hessian_ss = dx_ds.dot(dx_ds) + difference.dot(d2x_ds2);
                const Precision hessian_rs = dx_dr.dot(dx_ds) + difference.dot(d2x_drds);

                StaticMatrix<2, 2> hessian{};
                Vec2               gradient{};

                hessian << hessian_rr, hessian_rs,
                           hessian_rs, hessian_ss;

                gradient << gradient_r,
                            gradient_s;

                const Vec2 delta = hessian.ldlt().solve(-gradient);

                r += delta(0);
                s += delta(1);

                if (delta.norm() < tolerance) {
                    break;
                }
            }

            const Vec2 local{r, s};

            if (clip && !in_bounds(local)) {
                continue;
            }

            const Precision distance_squared = squared_distance(local);

            if (distance_squared < best_distance_squared) {
                best_local            = local;
                best_distance_squared = distance_squared;
            }
        }

        // Compare the interior projection against the closest boundary point
        if (clip) {
            const Vec2      boundary_local            = closest_point_on_boundary(global, coordinates);
            const Precision boundary_distance_squared = squared_distance(boundary_local);

            if (boundary_distance_squared < best_distance_squared) {
                best_local = boundary_local;
            }
        }

        return best_local;
    }

    /**
     * @brief Computes the normalized surface normal.
     */
    Vec3 normal(const Field& node_coords,
                const Vec2&  local) const override {
        const auto coordinates = node_coords_global(node_coords);
        const auto jac         = jacobian(coordinates, local(0), local(1));

        return jac.col(0).cross(jac.col(1)).normalized();
    }

    /**
     * @brief Integrates each shape function over the global surface.
     */
    DynamicVector shape_function_integral(const Field& node_coords) const override {
        const auto coordinates = node_coords_global(node_coords);

        std::function<StaticMatrix<N, 1>(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) {
                const auto shape = shape_function(r, s);
                const auto jac   = jacobian(coordinates, r, s);

                return shape * jac.col(0).cross(jac.col(1)).norm();
            };

        const auto integral = integration_scheme().integrate(integrand);
        return DynamicVector(integral);
    }

    /**
     * @brief Computes the global surface area.
     */
    Precision area(const Field& node_coords) const override {
        const auto coordinates = node_coords_global(node_coords);

        std::function<Precision(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) {
                const auto jac = jacobian(coordinates, r, s);
                return jac.col(0).cross(jac.col(1)).norm();
            };

        return integration_scheme().integrate(integrand);
    }

    /**
     * @brief Collects the global coordinates of the surface nodes.
     */
    StaticMatrix<N, 3> node_coords_global(const Field& node_coords) const {
        StaticMatrix<N, 3> coordinates{};

        for (Index local_id = 0; local_id < N; ++local_id) {
            const Index node_id = static_cast<Index>(nodeIds[local_id]);
            coordinates.row(local_id) = node_coords.row_vec3(node_id).transpose();
        }

        return coordinates;
    }

    /**
     * @brief Returns the global node identifiers.
     */
    const std::array<ID, N>& get_node_ids() const {
        return nodeIds;
    }

    /**
     * @brief Returns a mutable pointer to the global node identifiers.
     */
    ID* nodes() override {
        return nodeIds.data();
    }

    /**
     * @brief Applies a distributed global load to the surface nodes.
     */
    void apply_dload(const Field& node_coords,
                     Field&       node_loads,
                     Vec3         load) override {
        const auto nodal_contributions = shape_function_integral(node_coords);

        for (Index local_id = 0; local_id < n_nodes; ++local_id) {
            const ID node_id = nodeIds[local_id];

            for (Dim dof = 0; dof < 3; ++dof) {
                node_loads(node_id, dof) += nodal_contributions(local_id) * load(dof);
            }
        }
    }

    /**
     * @brief Applies a pressure load along the oriented surface normal.
     */
    void apply_pload(const Field& node_coords,
                     Field&       node_loads,
                     Precision    load) override {
        const auto coordinates = node_coords_global(node_coords);

        std::function<StaticMatrix<N, 3>(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) {
                const auto shape          = shape_function(r, s);
                const auto jac            = jacobian(coordinates, r, s);
                const Vec3 surface_vector = jac.col(0).cross(jac.col(1));

                return shape * (load * surface_vector).transpose();
            };

        const auto nodal_contributions = integration_scheme().integrate(integrand);

        for (Index local_id = 0; local_id < n_nodes; ++local_id) {
            const ID node_id = nodeIds[local_id];

            for (Dim dof = 0; dof < 3; ++dof) {
                node_loads(node_id, dof) += nodal_contributions(local_id, dof);
            }
        }
    }
};

} // namespace fem::model
