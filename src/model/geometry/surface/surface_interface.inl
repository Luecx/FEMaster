/**
 * @file surface_interface.inl
 * @brief Implements clipped triangle integration for SurfaceInterface.
 */

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace fem::model {
namespace surface_interface_detail {

/**
 * @brief Computes the signed two-dimensional cross product.
 */
inline Precision cross_2d(const Vec2& lhs, const Vec2& rhs) {
    return lhs(0) * rhs(1) - lhs(1) * rhs(0);
}

/**
 * @brief Appends the intersection of a polygon edge and a clipping boundary.
 */
inline void append_intersection(std::vector<Vec2>& clipped,
                                const Vec2&        previous,
                                const Vec2&        current,
                                Precision          previous_side,
                                Precision          current_side,
                                Precision          tolerance) {
    const Precision denominator = previous_side - current_side;
    if (std::abs(denominator) <= tolerance) {
        return;
    }

    const Precision alpha = previous_side / denominator;
    clipped.push_back(previous + alpha * (current - previous));
}

/**
 * @brief Clips a polygon against the convex local surface domain.
 */
inline std::vector<Vec2> clip_polygon_to_domain(std::vector<Vec2>       polygon,
                                                const std::vector<Vec2>& domain,
                                                Precision                tolerance) {
    for (std::size_t edge_id = 0; edge_id < domain.size(); ++edge_id) {
        if (polygon.empty()) {
            return {};
        }

        const Vec2 edge_start     = domain[edge_id];
        const Vec2 edge_end       = domain[(edge_id + 1) % domain.size()];
        const Vec2 edge_direction = edge_end - edge_start;

        std::vector<Vec2> clipped{};
        clipped.reserve(polygon.size() + 1);

        Vec2      previous        = polygon.back();
        Precision previous_side   = cross_2d(edge_direction, previous - edge_start);
        bool      previous_inside = previous_side >= -tolerance;

        for (const Vec2& current : polygon) {
            const Precision current_side   = cross_2d(edge_direction, current - edge_start);
            const bool      current_inside = current_side >= -tolerance;

            if (previous_inside && current_inside) {
                clipped.push_back(current);
            } else if (previous_inside && !current_inside) {
                append_intersection(
                    clipped,
                    previous,
                    current,
                    previous_side,
                    current_side,
                    tolerance
                );
            } else if (!previous_inside && current_inside) {
                append_intersection(
                    clipped,
                    previous,
                    current,
                    previous_side,
                    current_side,
                    tolerance
                );

                clipped.push_back(current);
            }

            previous        = current;
            previous_side   = current_side;
            previous_inside = current_inside;
        }

        polygon = std::move(clipped);
    }

    return polygon;
}

/**
 * @brief Removes consecutive and closing duplicate polygon vertices.
 */
inline std::vector<Vec2> clean_polygon(const std::vector<Vec2>& polygon,
                                       Precision                tolerance) {
    std::vector<Vec2> cleaned{};
    cleaned.reserve(polygon.size());

    for (const Vec2& point : polygon) {
        if (cleaned.empty() || (point - cleaned.back()).norm() > tolerance) {
            cleaned.push_back(point);
        }
    }

    if (cleaned.size() > 1 &&
        (cleaned.front() - cleaned.back()).norm() <= tolerance) {
        cleaned.pop_back();
    }

    return cleaned;
}

} // namespace surface_interface_detail

template<typename F>
void SurfaceInterface::integrate_local_triangle(const std::array<Vec2, 3>&    triangle,
                                                const Field&                  node_coords,
                                                const quadrature::Quadrature& scheme,
                                                const F&                      integrand) const {
    constexpr Precision tolerance = Precision(1e-12);

    // Clip the triangle to the local surface domain
    std::vector<Vec2> polygon = {
        triangle[0],
        triangle[1],
        triangle[2]
    };

    polygon = surface_interface_detail::clip_polygon_to_domain(
        std::move(polygon),
        local_domain_polygon(),
        tolerance
    );

    // Remove duplicate vertices introduced by clipping
    polygon = surface_interface_detail::clean_polygon(polygon, tolerance);

    if (polygon.size() < 3) {
        return;
    }

    // Triangulate and integrate the clipped polygon
    const Vec2& origin = polygon.front();

    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        const Vec2 edge_r = polygon[i]     - origin;
        const Vec2 edge_s = polygon[i + 1] - origin;

        const Precision triangle_jacobian =
            std::abs(surface_interface_detail::cross_2d(edge_r, edge_s));

        if (triangle_jacobian <= tolerance) {
            continue;
        }

        for (Index q = 0; q < scheme.count(); ++q) {
            const auto point = scheme.get_point(q);

            const Vec2 local  = origin + point.r * edge_r + point.s * edge_s;
            const Vec3 global = local_to_global(local, node_coords);

            const Precision weight = point.w
                                   * triangle_jacobian
                                   * jacobian_measure(node_coords, local);

            integrand(local, global, weight);
        }
    }
}

} // namespace fem::model