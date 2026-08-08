/**
 * @file nagata_surface_project.cpp
 * @brief Implements global and tracked closest-point projection on Nagata surfaces.
 *
 * Global projection uses the nearest reconstructed vertex as a seed source.
 * Tracked projection solves the closest-point equations with a consistent
 * Newton Hessian and transfers barycentric edge locations across adjacent
 * patches. True surface boundaries use a one-dimensional edge projection.
 *
 * @see NagataSurface
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "nagata_surface.h"

#include "../../../core/logging.h"

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem::model {

/**
 * @brief Projects a finite global point using all patches at its nearest vertex.
 *
 * The vertex BVH first identifies the nearest reconstructed node. Every incident
 * patch is then seeded slightly inside its reference triangle so evaluation does
 * not begin at the removable rational singularity of a patch vertex. Tracked
 * projection may walk across further patches; the candidate with minimum
 * squared geometric distance is returned.
 *
 * @param point Finite point in global Cartesian coordinates.
 * @return Closest reconstructed location found from the incident patch fan.
 */
NagataSurface::Location NagataSurface::project(const Vec3& point) const {
    logging::error(point.allFinite(),
        "NagataSurface::project requires a finite point");
    logging::error(vertex_bvh_.valid(),
        "NagataSurface::project requires a valid vertex BVH");

    const ID nearest_vertex_id = vertex_bvh_.query_nearest(point);

    logging::error(nearest_vertex_id >= 0,
        "NagataSurface::project could not determine a nearest surface vertex");

    const Index nearest_vertex = static_cast<Index>(nearest_vertex_id);

    logging::error(nearest_vertex < static_cast<Index>(vertices_.size()),
        "NagataSurface::project BVH returned an invalid vertex ID");

    const Vertex& vertex = vertices_[static_cast<std::size_t>(nearest_vertex)];

    logging::error(!vertex.patches.empty(),
        "NagataSurface::project nearest vertex has no incident patches");

    constexpr Precision seed_offset = Precision(1e-6);

    Precision best_distance = std::numeric_limits<Precision>::max();
    Location  best_location;

    for (const nagata::PatchID patch_id : vertex.patches) {
        const Patch& patch = patches_[static_cast<std::size_t>(patch_id)];
        Index local_vertex = 3;

        for (Index i = 0; i < 3; ++i) {
            if (patch.vertices[static_cast<std::size_t>(i)] == nearest_vertex) {
                local_vertex = i;
                break;
            }
        }

        logging::error(local_vertex < 3,
            "NagataSurface::project vertex incidence is inconsistent");

        Location initial_guess;
        initial_guess.patch = patch_id;

        // Reduced coordinates equal beta_1 and beta_2. Retain the selected
        // vertex as dominant while moving a small distance into the patch.
        if (local_vertex == 0) {
            initial_guess.local << seed_offset, seed_offset;
        } else if (local_vertex == 1) {
            initial_guess.local <<
                Precision(1) - Precision(2)*seed_offset,
                seed_offset;
        } else {
            initial_guess.local <<
                seed_offset,
                Precision(1) - Precision(2)*seed_offset;
        }

        const Location candidate = project(point, initial_guess);
        const Evaluation state   = evaluate(candidate);
        const Precision distance = (state.position - point).squaredNorm();

        if (distance < best_distance) {
            best_distance = distance;
            best_location = candidate;
        }
    }

    logging::error(best_location.patch != invalid_patch,
        "NagataSurface::project failed to produce a surface location");

    return best_location;
}

/**
 * @brief Continues closest-point projection from an existing surface location.
 *
 * Closest-point coordinates satisfy `x_,alpha dot (x - p) = 0`. The Newton
 * Hessian combines the surface metric with the residual contraction of the
 * second derivatives. A singular consistent Hessian falls back to the positive
 * semidefinite Gauss-Newton metric.
 *
 * Steps leaving a patch are clipped at the first crossed barycentric edge. An
 * internal edge transfers the same shared-vertex barycentric weights to the
 * neighbor. A boundary edge instead solves its one-dimensional closest-point
 * equation with a clamped Newton iteration.
 *
 * @param point Finite global point to project.
 * @param initial_guess Valid location used as the Newton seed.
 * @return Projected location after local patch walking or boundary projection.
 */
NagataSurface::Location NagataSurface::project(
    const Vec3&     point,
    const Location& initial_guess) const {
    logging::error(point.allFinite(),
        "NagataSurface::project requires a finite point");
    logging::error(initial_guess.patch < static_cast<nagata::PatchID>(patches_.size()),
        "NagataSurface::project received an invalid initial patch");
    logging::error(initial_guess.local.allFinite(),
        "NagataSurface::project received invalid initial coordinates");

    constexpr Index maximum_iterations          = 30;
    constexpr Index maximum_boundary_iterations = 20;

    const Precision tolerance =
        std::sqrt(std::numeric_limits<Precision>::epsilon());

    Location location = initial_guess;

    for (Index iteration = 0; iteration < maximum_iterations; ++iteration) {
        // -------------------------------------------------------------
        // Evaluate the closest-point residual and Newton gradient
        // -------------------------------------------------------------

        const Evaluation state = evaluate(location);
        const Vec3 residual     = state.position - point;
        const Vec3 tangent_xi   = state.jacobian.col(0);
        const Vec3 tangent_eta  = state.jacobian.col(1);

        Vec2 gradient;
        gradient(0) = tangent_xi.dot(residual);
        gradient(1) = tangent_eta.dot(residual);

        const Precision convergence_scale =
            Precision(1) + state.jacobian.norm() * residual.norm();

        if (gradient.norm() <= tolerance * convergence_scale) {
            return location;
        }

        // -------------------------------------------------------------
        // Assemble and solve the consistent closest-point Hessian
        // -------------------------------------------------------------

        StaticMatrix<2, 2> hessian;
        hessian(0, 0) = tangent_xi.dot(tangent_xi)
                      + residual.dot(state.d2_xixi);
        hessian(0, 1) = tangent_xi.dot(tangent_eta)
                      + residual.dot(state.d2_xieta);
        hessian(1, 0) = hessian(0, 1);
        hessian(1, 1) = tangent_eta.dot(tangent_eta)
                      + residual.dot(state.d2_etaeta);

        Vec2 increment;
        Eigen::FullPivLU<StaticMatrix<2, 2>> hessian_lu(hessian);

        if (hessian_lu.isInvertible()) {
            increment = -hessian_lu.solve(gradient);
        } else {
            const StaticMatrix<2, 2> metric =
                state.jacobian.transpose() * state.jacobian;
            Eigen::FullPivLU<StaticMatrix<2, 2>> metric_lu(metric);

            logging::error(metric_lu.isInvertible(),
                "NagataSurface::project encountered a degenerate surface metric");

            increment = -metric_lu.solve(gradient);
        }

        logging::error(increment.allFinite(),
            "NagataSurface::project produced an invalid Newton increment");

        const Vec2 candidate_local = location.local + increment;

        // Express both ends of the Newton segment in full barycentric form.
        Vec3 beta_current;
        beta_current <<
            Precision(1) - location.local(0) - location.local(1),
            location.local(0),
            location.local(1);

        Vec3 beta_candidate;
        beta_candidate <<
            Precision(1) - candidate_local(0) - candidate_local(1),
            candidate_local(0),
            candidate_local(1);

        if (beta_candidate.minCoeff() >= Precision(0)) {
            location.local = candidate_local;

            if (increment.norm() <= tolerance) {
                return location;
            }

            continue;
        }

        // -------------------------------------------------------------
        // Locate the first edge crossed by the Newton segment
        // -------------------------------------------------------------

        Precision crossing = Precision(1);
        Index     edge     = 3;

        for (Index i = 0; i < 3; ++i) {
            if (beta_candidate(static_cast<Eigen::Index>(i)) >= Precision(0)) {
                continue;
            }

            const Precision current = beta_current(static_cast<Eigen::Index>(i));
            const Precision candidate = beta_candidate(static_cast<Eigen::Index>(i));
            const Precision denominator = current - candidate;

            if (denominator <= Precision(0)) {
                continue;
            }

            const Precision alpha = current / denominator;

            if (alpha < crossing) {
                crossing = alpha;
                edge     = i;
            }
        }

        logging::error(edge < 3,
            "NagataSurface::project could not identify a crossed patch edge");

        Vec3 beta_edge =
            beta_current + crossing * (beta_candidate - beta_current);
        beta_edge(static_cast<Eigen::Index>(edge)) = Precision(0);

        const Patch& current_patch = patches_[static_cast<std::size_t>(location.patch)];
        const Neighbor& neighbor   = current_patch.neighbors[static_cast<std::size_t>(edge)];

        // -------------------------------------------------------------
        // Solve a one-dimensional projection on a true boundary edge
        // -------------------------------------------------------------

        if (neighbor.patch == invalid_patch) {
            const Index j = (edge + 1) % 3;
            const Index k = (edge + 2) % 3;

            const std::array<Vec2, 3> local_vertices{
                Vec2(Precision(0), Precision(0)),
                Vec2(Precision(1), Precision(0)),
                Vec2(Precision(0), Precision(1))
            };

            const Vec2 local_j  = local_vertices[static_cast<std::size_t>(j)];
            const Vec2 local_k  = local_vertices[static_cast<std::size_t>(k)];
            const Vec2 direction = local_k - local_j;

            const Precision edge_sum =
                beta_edge(static_cast<Eigen::Index>(j))
                + beta_edge(static_cast<Eigen::Index>(k));

            Precision edge_parameter = edge_sum > Precision(0)
                ? beta_edge(static_cast<Eigen::Index>(k)) / edge_sum
                : Precision(0.5);

            edge_parameter = std::clamp(edge_parameter, Precision(0), Precision(1));

            Location boundary_location;
            boundary_location.patch = location.patch;

            for (Index boundary_iteration = 0;
                 boundary_iteration < maximum_boundary_iterations;
                 ++boundary_iteration) {
                boundary_location.local =
                    (Precision(1) - edge_parameter) * local_j
                    + edge_parameter * local_k;

                const Evaluation boundary_state = evaluate(boundary_location);
                const Vec3 boundary_residual = boundary_state.position - point;
                const Vec3 tangent = boundary_state.jacobian * direction;
                const Vec3 second_derivative =
                      direction(0)*direction(0) * boundary_state.d2_xixi
                    + Precision(2)*direction(0)*direction(1) * boundary_state.d2_xieta
                    + direction(1)*direction(1) * boundary_state.d2_etaeta;

                const Precision gradient_1d = tangent.dot(boundary_residual);
                const Precision hessian_1d = tangent.dot(tangent)
                                             + boundary_residual.dot(second_derivative);
                const Precision boundary_scale =
                    Precision(1) + tangent.norm() * boundary_residual.norm();

                if (std::abs(gradient_1d) <= tolerance * boundary_scale) {
                    return boundary_location;
                }

                Precision effective_hessian = hessian_1d;

                if (std::abs(effective_hessian)
                    <= std::numeric_limits<Precision>::epsilon()) {
                    effective_hessian = tangent.squaredNorm();
                }

                if (std::abs(effective_hessian)
                    <= std::numeric_limits<Precision>::epsilon()) {
                    return boundary_location;
                }

                const Precision next_parameter = std::clamp(
                    edge_parameter - gradient_1d / effective_hessian,
                    Precision(0),
                    Precision(1));

                if (std::abs(next_parameter - edge_parameter) <= tolerance) {
                    edge_parameter = next_parameter;
                    boundary_location.local =
                        (Precision(1) - edge_parameter) * local_j
                        + edge_parameter * local_k;

                    return boundary_location;
                }

                edge_parameter = next_parameter;
            }

            logging::warning(false,
                "NagataSurface::project boundary projection reached the iteration limit");

            return boundary_location;
        }

        // -------------------------------------------------------------
        // Transfer shared-edge weights to the neighboring smooth patch
        // -------------------------------------------------------------

        logging::error(neighbor.patch < static_cast<nagata::PatchID>(patches_.size()),
            "NagataSurface::project contains invalid patch adjacency");

        const Patch& next_patch = patches_[static_cast<std::size_t>(neighbor.patch)];
        Vec3 next_beta = Vec3::Zero();

        const Index source_j = (edge + 1) % 3;
        const Index source_k = (edge + 2) % 3;
        Index matched_vertices = 0;

        for (const Index source_local : {source_j, source_k}) {
            const Index vertex =
                current_patch.vertices[static_cast<std::size_t>(source_local)];
            bool matched = false;

            for (Index target_local = 0; target_local < 3; ++target_local) {
                if (next_patch.vertices[static_cast<std::size_t>(target_local)] != vertex) {
                    continue;
                }

                next_beta(static_cast<Eigen::Index>(target_local)) =
                    beta_edge(static_cast<Eigen::Index>(source_local));
                matched = true;
                ++matched_vertices;
                break;
            }

            logging::error(matched,
                "NagataSurface::project neighboring patches do not share the expected edge");
        }

        logging::error(matched_vertices == 2,
            "NagataSurface::project failed to transfer patch-edge coordinates");

        location.patch = neighbor.patch;
        location.local << next_beta(1), next_beta(2);
    }

    logging::warning(false,
        "NagataSurface::project reached the maximum number of Newton iterations");

    return location;
}

} // namespace fem::model
