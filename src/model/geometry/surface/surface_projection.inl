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
    constexpr Index     max_iterations             = 40;
    constexpr Index     max_line_search_iterations = 20;
    constexpr Precision step_tolerance             = Precision(1e-12);
    constexpr Precision gradient_tolerance         = Precision(1e-12);
    constexpr Precision armijo_constant            = Precision(1e-4);
    constexpr Precision minimum_step_length        = Precision(1e-8);
    constexpr Precision metric_regularization      = Precision(1e-14);

    const auto coordinates = node_coords_global(node_coords);

    const auto objective = [&](const Vec2& local) {
        const Vec3 position = interpolate(coordinates, local(0), local(1));
        return Precision(0.5) * (position - global).squaredNorm();
    };

    Vec2      best_local     = Vec2::Zero();
    Precision best_objective = std::numeric_limits<Precision>::max();

    const auto project_initial_guess = [&](const Vec2& initial_guess) {
        Vec2      local             = initial_guess;
        Precision current_objective = objective(local);

        if (!std::isfinite(current_objective)) {
            return;
        }

        for (Index iteration = 0; iteration < max_iterations; ++iteration) {
            const Precision r = local(0);
            const Precision s = local(1);

            const auto shape             = shape_function(r, s);
            const auto derivative        = shape_derivative(r, s);
            const auto second_derivative = shape_second_derivative(r, s);

            const Vec3               position = coordinates.transpose() * shape;
            const StaticMatrix<3, 2> first    = coordinates.transpose() * derivative;
            const StaticMatrix<3, 3> second   = coordinates.transpose() * second_derivative;

            const Vec3 dx_dr = first.col(0);
            const Vec3 dx_ds = first.col(1);

            const Vec3 d2x_dr2  = second.col(0);
            const Vec3 d2x_ds2  = second.col(1);
            const Vec3 d2x_drds = second.col(2);

            const Vec3 difference = position - global;

            Vec2 gradient;
            gradient << difference.dot(dx_dr),
                        difference.dot(dx_ds);

            if (!gradient.allFinite() || gradient.norm() < gradient_tolerance) {
                break;
            }

            StaticMatrix<2, 2> hessian;
            hessian << dx_dr.dot(dx_dr) + difference.dot(d2x_dr2),
                       dx_dr.dot(dx_ds) + difference.dot(d2x_drds),
                       dx_dr.dot(dx_ds) + difference.dot(d2x_drds),
                       dx_ds.dot(dx_ds) + difference.dot(d2x_ds2);

            Vec2 direction = Vec2::Zero();

            if (hessian.allFinite()) {
                const Eigen::LDLT<StaticMatrix<2, 2>> decomposition(hessian);

                if (decomposition.info() == Eigen::Success) {
                    direction = decomposition.solve(-gradient);
                }
            }

            if (!direction.allFinite() || gradient.dot(direction) >= Precision(0)) {
                StaticMatrix<2, 2> metric = first.transpose() * first;
                metric.diagonal().array() += metric_regularization;

                const Eigen::LDLT<StaticMatrix<2, 2>> decomposition(metric);

                if (decomposition.info() == Eigen::Success) {
                    direction = decomposition.solve(-gradient);
                }
            }

            if (!direction.allFinite() || gradient.dot(direction) >= Precision(0)) {
                direction = -gradient;
            }

            const Precision directional_derivative = gradient.dot(direction);

            if (!std::isfinite(directional_derivative) ||
                directional_derivative >= Precision(0)) {
                break;
            }

            Precision step_length       = Precision(1);
            Vec2      accepted_local     = local;
            Precision accepted_objective = current_objective;
            bool      accepted           = false;

            for (Index line_search_iteration = 0;
                 line_search_iteration < max_line_search_iterations;
                 ++line_search_iteration) {
                const Vec2 trial_local = local + step_length * direction;

                if (clip && !in_bounds(trial_local)) {
                    step_length *= Precision(0.5);
                    continue;
                }

                const Precision trial_objective = objective(trial_local);
                const Precision armijo_bound =
                    current_objective +
                    armijo_constant * step_length * directional_derivative;

                if (std::isfinite(trial_objective) &&
                    trial_objective <= armijo_bound) {
                    accepted_local     = trial_local;
                    accepted_objective = trial_objective;
                    accepted           = true;
                    break;
                }

                step_length *= Precision(0.5);
            }

            if (!accepted || step_length < minimum_step_length) {
                break;
            }

            const Vec2 applied_step = accepted_local - local;

            local             = accepted_local;
            current_objective = accepted_objective;

            if (applied_step.norm() < step_tolerance) {
                break;
            }
        }

        if (clip && !in_bounds(local)) {
            return;
        }

        const Precision final_objective = objective(local);

        if (std::isfinite(final_objective) &&
            final_objective < best_objective) {
            best_local     = local;
            best_objective = final_objective;
        }
    };

    if constexpr (N == 3) {
        const std::array initial_guesses{
            Vec2{Precision(1.0 / 3.0), Precision(1.0 / 3.0)},
            Vec2{Precision(0.10), Precision(0.10)},
            Vec2{Precision(0.80), Precision(0.10)},
            Vec2{Precision(0.10), Precision(0.80)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else if constexpr (N == 4) {
        const std::array initial_guesses{
            Vec2{Precision( 0.00), Precision( 0.00)},
            Vec2{Precision(-0.50), Precision(-0.50)},
            Vec2{Precision( 0.50), Precision(-0.50)},
            Vec2{Precision( 0.50), Precision( 0.50)},
            Vec2{Precision(-0.50), Precision( 0.50)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else if constexpr (N == 6) {
        const std::array initial_guesses{
            Vec2{Precision(1.0 / 3.0), Precision(1.0 / 3.0)},
            Vec2{Precision(0.10), Precision(0.10)},
            Vec2{Precision(0.80), Precision(0.10)},
            Vec2{Precision(0.10), Precision(0.80)},
            Vec2{Precision(0.45), Precision(0.10)},
            Vec2{Precision(0.45), Precision(0.45)},
            Vec2{Precision(0.10), Precision(0.45)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    } else {
        const std::array initial_guesses{
            Vec2{Precision( 0.00), Precision( 0.00)},
            Vec2{Precision(-0.50), Precision(-0.50)},
            Vec2{Precision( 0.50), Precision(-0.50)},
            Vec2{Precision( 0.50), Precision( 0.50)},
            Vec2{Precision(-0.50), Precision( 0.50)},
            Vec2{Precision( 0.00), Precision(-0.80)},
            Vec2{Precision( 0.80), Precision( 0.00)},
            Vec2{Precision( 0.00), Precision( 0.80)},
            Vec2{Precision(-0.80), Precision( 0.00)}
        };

        for (const Vec2& initial_guess : initial_guesses) {
            project_initial_guess(initial_guess);
        }
    }

    if (clip) {
        const Vec2      boundary_local     = closest_point_on_boundary(global, coordinates);
        const Precision boundary_objective = objective(boundary_local);

        if (std::isfinite(boundary_objective) &&
            boundary_objective < best_objective) {
            best_local     = boundary_local;
            best_objective = boundary_objective;
        }
    }

    return best_local;
}
} // namespace fem::model
