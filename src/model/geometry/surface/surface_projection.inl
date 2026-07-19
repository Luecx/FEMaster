/**
 * @file surface_projection.inl
 * @brief Implements global-to-local projection for finite-element surfaces.
 *
 * @see surface.h
 */

#pragma once

#include <Eigen/Cholesky>

#include <array>
#include <limits>

namespace fem::model {

/**
 * @brief Projects a global point onto the finite-element surface.
 *
 * The method minimizes the squared physical distance between the requested
 * global point and the isoparametric surface. Newton iterations are started
 * from several topology-dependent natural coordinates to reduce the risk of
 * converging to an unsuitable stationary point.
 *
 * When `clip` is enabled, unconstrained solutions outside the natural element
 * domain are rejected. The best valid interior solution is then compared with
 * the closest point on the element boundary.
 *
 * @param global Global point to project.
 * @param node_coords Global nodal coordinate field.
 * @param clip Whether the returned point must lie inside the natural domain.
 *
 * @return Natural coordinates of the closest projection found.
 */
template<Index N>
Vec2 Surface<N>::global_to_local(
    const Vec3&  global,
    const Field& node_coords,
    bool         clip
) const {
    constexpr Index     max_iterations = 32;
    constexpr Precision tolerance      = Precision(1e-12);

    // collect the fixed-size coordinate matrix once because it is reused by
    // every initial guess and every Newton iteration
    const auto coordinates = node_coords_global(node_coords);

    // evaluate the physical squared distance associated with natural
    // coordinates. squared distances are sufficient for all comparisons.
    const auto squared_distance = [&](const Vec2& local) {
        const Vec3 position = interpolate(coordinates, local(0), local(1));
        return (position - global).squaredNorm();
    };

    Vec2      best_local            = Vec2::Zero();
    Precision best_distance_squared = std::numeric_limits<Precision>::max();

    // project one initial guess with Newton's method and retain it when it
    // represents the best valid stationary point found so far
    const auto project_initial_guess = [&](const Vec2& initial_guess) {
        Precision r = initial_guess(0);
        Precision s = initial_guess(1);

        for (Index iteration = 0; iteration < max_iterations; ++iteration) {
            // evaluate the interpolation and its first and second derivatives
            // at the current natural coordinates
            const auto shape             = shape_function(r, s);
            const auto derivative        = shape_derivative(r, s);
            const auto second_derivative = shape_second_derivative(r, s);

            // map the interpolation quantities into physical space
            const Vec3               position = coordinates.transpose() * shape;
            const StaticMatrix<3, 2> first    = coordinates.transpose() * derivative;
            const StaticMatrix<3, 3> second   = coordinates.transpose() * second_derivative;

            // the first derivative columns are the two physical tangent vectors
            const Vec3 dx_dr = first.col(0);
            const Vec3 dx_ds = first.col(1);

            // the second derivative columns follow the common ordering
            // d²x/dr², d²x/ds², and d²x/(dr ds)
            const Vec3 d2x_dr2  = second.col(0);
            const Vec3 d2x_ds2  = second.col(1);
            const Vec3 d2x_drds = second.col(2);

            const Vec3 difference = position - global;

            // gradient of one half of the squared-distance objective
            const Precision gradient_r = difference.dot(dx_dr);
            const Precision gradient_s = difference.dot(dx_ds);

            // exact Hessian of one half of the squared-distance objective
            const Precision hessian_rr = dx_dr.dot(dx_dr) + difference.dot(d2x_dr2);
            const Precision hessian_ss = dx_ds.dot(dx_ds) + difference.dot(d2x_ds2);
            const Precision hessian_rs = dx_dr.dot(dx_ds) + difference.dot(d2x_drds);

            StaticMatrix<2, 2> hessian;
            Vec2               gradient;

            hessian << hessian_rr, hessian_rs,
                       hessian_rs, hessian_ss;

            gradient << gradient_r,
                        gradient_s;

            // solve the two-dimensional Newton system for the natural-coordinate
            // correction
            const Vec2 delta = hessian.ldlt().solve(-gradient);

            r += delta(0);
            s += delta(1);

            // the natural-coordinate correction directly provides the local
            // convergence measure
            if (delta.norm() < tolerance) {
                break;
            }
        }

        const Vec2 local{r, s};

        // when clipping is requested, stationary points outside the finite
        // natural element domain must not be considered
        if (clip && !in_bounds(local)) {
            return;
        }

        const Precision distance_squared = squared_distance(local);

        // retain the closest valid projection found across all initial guesses
        if (distance_squared < best_distance_squared) {
            best_local            = local;
            best_distance_squared = distance_squared;
        }
    };

    // use topology-specific fixed-size initial guesses without dynamic memory
    // allocation
    if constexpr (N == 3) {
        const std::array initial_guesses{
            Vec2{Precision(0.25), Precision(0.25)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else if constexpr (N == 4) {
        const std::array initial_guesses{
            Vec2{Precision(-1), Precision(-1)},
            Vec2{Precision( 1), Precision(-1)},
            Vec2{Precision( 1), Precision( 1)},
            Vec2{Precision(-1), Precision( 1)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else if constexpr (N == 6) {
        const std::array initial_guesses{
            Vec2{Precision(0), Precision(0)},
            Vec2{Precision(1), Precision(0)},
            Vec2{Precision(0), Precision(1)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else {
        const std::array initial_guesses{
            Vec2{Precision(-1), Precision(-1)},
            Vec2{Precision( 1), Precision(-1)},
            Vec2{Precision( 1), Precision( 1)},
            Vec2{Precision(-1), Precision( 1)},
            Vec2{Precision( 0), Precision(-1)},
            Vec2{Precision( 1), Precision( 0)},
            Vec2{Precision( 0), Precision( 1)},
            Vec2{Precision(-1), Precision( 0)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    }

    // a constrained projection can lie either in the element interior or on
    // its boundary, so compare the best interior candidate with the closest
    // boundary projection
    if (clip) {
        const Vec2 boundary_local = closest_point_on_boundary(global, coordinates);
        const Precision boundary_distance_squared = squared_distance(boundary_local);

        if (boundary_distance_squared < best_distance_squared) {
            best_local = boundary_local;
        }
    }

    return best_local;
}

} // namespace fem::model