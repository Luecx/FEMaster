/**
 * @file surface_polygon.inl
 * @brief Implements the fixed-capacity polygon for surface operations in 2D.
 *
 * The implementation stores vertices in fixed-size arrays and evaluates
 * orientation through the signed shoelace area. Convex intersections are
 * computed by Sutherland-Hodgman clipping with alternating fixed-capacity
 * buffers and a small geometric tolerance.
 *
 * @see SurfacePolygon
 * @see surface_polygon.h
 */

#pragma once

#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"

#include <cmath>
#include <utility>

namespace fem {

/**
 * Constructs a polygon from an initializer list of vertices.
 *
 * The supplied vertices are copied in their given order. The list must not
 * contain more vertices than the compile-time capacity `MP`.
 *
 * @param points Vertices of the polygon in boundary order.
 */
template<std::size_t MP>
SurfacePolygon<MP>::SurfacePolygon(std::initializer_list<Vec2> points) {
    // The number of supplied points must fit into the fixed-capacity storage
    logging::error(points.size() <= MP,
        "SurfacePolygon initializer contains more points than its fixed capacity");

    // Copy the supplied points into the internal array and update the actual polygon size
    for (const Vec2& point : points) {
        points_[size_++] = point;
    }
}

/**
 * Returns the number of active vertices in the polygon.
 *
 * The fixed-capacity backing array may contain unused entries; only the first
 * `size()` entries belong to the polygon.
 *
 * @return Number of active vertices.
 */
template<std::size_t MP>
std::size_t SurfacePolygon<MP>::size() const {
    return size_;
}

/**
 * Returns a checked read-only vertex reference.
 *
 * @param index Zero-based active vertex index.
 * @return Read-only reference to the requested vertex.
 */
template<std::size_t MP>
const Vec2& SurfacePolygon<MP>::operator[](std::size_t index) const {
    // Only the first size_ entries of points_ belong to the actual polygon
    logging::error(index < size_,
        "SurfacePolygon point index is out of bounds");

    return points_[index];
}

/**
 * Returns a checked mutable vertex reference.
 *
 * @param index Zero-based active vertex index.
 * @return Mutable reference to the requested vertex.
 */
template<std::size_t MP>
Vec2& SurfacePolygon<MP>::operator[](std::size_t index) {
    // Only the first size_ entries of points_ belong to the actual polygon
    logging::error(index < size_,
        "SurfacePolygon point index is out of bounds");

    return points_[index];
}

/**
 * Appends one vertex to the fixed-capacity polygon.
 *
 * The vertex becomes the final active vertex and the polygon boundary order
 * is preserved.
 *
 * @param point Vertex to append.
 */
template<std::size_t MP>
void SurfacePolygon<MP>::push_back(const Vec2& point) {
    // SurfacePolygon has fixed storage and must never allocate dynamically
    logging::error(size_ < MP,
        "SurfacePolygon fixed capacity exceeded while adding a point");

    points_[size_++] = point;
}

/**
 * Reverses the polygon orientation while retaining its first vertex.
 *
 * Keeping the first vertex fixed makes the operation deterministic while
 * changing clockwise ordering to counter-clockwise ordering and vice versa.
 * Polygons with fewer than three vertices have no meaningful orientation.
 */
template<std::size_t MP>
void SurfacePolygon<MP>::flip() {
    // Polygons containing fewer than three points do not have a meaningful orientation
    if (size_ < 3) {
        return;
    }

    // Reverse the winding while retaining the first point. keeping the first
    // Point fixed makes the operation deterministic with respect to the polygon
    // Anchor while still changing between clockwise and counter-clockwise order.
    for (std::size_t left = 1, right = size_ - 1; left < right; ++left, --right) {
        std::swap(points_[left], points_[right]);
    }
}

/**
 * Evaluates twice the signed polygon area using the shoelace formula.
 *
 * The result is positive for counter-clockwise ordering and negative for
 * clockwise ordering. Degenerate polygons naturally produce zero.
 *
 * @return Twice the signed polygon area.
 */
template<std::size_t MP>
Precision SurfacePolygon<MP>::area2() const {
    Precision result = Precision(0);

    // Evaluate the shoelace formula over all polygon edges. the modulo closes
    // The polygon by connecting the final point back to the first point.
    for (std::size_t i = 0; i < size_; ++i) {
        const Vec2& a = points_[i];
        const Vec2& b = points_[(i + 1) % size_];

        result += a(0) * b(1) - a(1) * b(0);
    }

    // The sign contains the polygon orientation:
    // Positive for counter-clockwise and negative for clockwise ordering
    return result;
}

/**
 * Computes the arithmetic mean of all active polygon vertices.
 *
 * This is the vertex average, not an area-weighted centroid. The polygon must
 * contain at least one vertex.
 *
 * @return Arithmetic mean of the active vertices.
 */
template<std::size_t MP>
Vec2 SurfacePolygon<MP>::mid() const {
    // An empty polygon does not provide a meaningful average point
    logging::error(size_ > 0,
        "SurfacePolygon cannot compute the midpoint of an empty polygon");
    
    Vec2 result{0,0};

    // Just sum over all points and then divide by amount of points
    for (std::size_t i = 0; i < size_; ++i) {
        result += points_[i];
    }

    // Divide to adjust for summation above
    return result / size_;
}

/**
 * Returns the non-negative geometric area of the polygon.
 *
 * The signed shoelace result is converted to an unsigned area by taking its
 * absolute value.
 *
 * @return Non-negative polygon area.
 */
template<std::size_t MP>
Precision SurfacePolygon<MP>::area() const {
    // Area2 is signed, while the physical polygon area is always non-negative
    return Precision(0.5) * std::abs(area2());
}

/**
 * Tests whether the active vertices are ordered counter-clockwise.
 *
 * Empty and degenerate polygons are treated as counter-clockwise because their
 * signed area is non-negative.
 *
 * @return `true` when the signed area is non-negative.
 */
template<std::size_t MP>
bool SurfacePolygon<MP>::is_ccw() const {
    // A non-negative signed area indicates counter-clockwise vertex ordering
    return area2() >= Precision(0);
}

/**
 * Computes the intersection of two convex polygons.
 *
 * Both polygons are normalized to counter-clockwise ordering and the subject
 * polygon is successively clipped against every edge of `clip_polygon` using
 * the Sutherland-Hodgman algorithm. Two fixed-capacity buffers are alternated
 * between clipping passes, so the operation performs no dynamic allocation.
 * A small tolerance classifies points near a clipping edge as inside and
 * suppresses unstable intersections of nearly parallel edges.
 *
 * @tparam MPO Capacity of the clipping polygon.
 * @param clip_polygon Convex polygon defining the clipping region.
 * @return Fixed-capacity polygon containing the convex intersection.
 */
template<std::size_t MP>
template<std::size_t MPO>
SurfacePolygon<MP + MPO> SurfacePolygon<MP>::intersection(
    const SurfacePolygon<MPO>& clip_polygon
) const {
    // The intersection of two convex polygons with MP and MPO vertices can
    // Contain at most MP + MPO vertices
    constexpr std::size_t capacity  = MP + MPO;
    constexpr Precision   tolerance = Precision(1e-12);

    // A valid polygonal intersection requires both inputs to contain at least
    // Three points
    if (size_ < 3 || clip_polygon.size() < 3) {
        return {};
    }

    // Sutherland-Hodgman interprets the left side of every clipping edge as
    // The interior half-plane. therefore the clipping polygon has to be ccw.
    // We normalize both polygons so the resulting intersection is also ccw.
    SurfacePolygon<MP>  subject = *this;
    SurfacePolygon<MPO> clip    = clip_polygon;

    if (!subject.is_ccw()) {
        subject.flip();
    }

    if (!clip.is_ccw()) {
        clip.flip();
    }

    // Clipping is performed repeatedly. two fixed-capacity buffers allow us
    // To alternate between input and output without dynamic allocation.
    std::array<Vec2, capacity> buffer_a{};
    std::array<Vec2, capacity> buffer_b{};

    // The first clipping pass reads directly from the normalized subject
    // Polygon. later passes read from the temporary clipping buffers.
    const Vec2* polygon      = subject.points_.data();
    Vec2*       clipped      = buffer_a.data();
    std::size_t polygon_size = subject.size_;

    // Signed two-dimensional cross product. for a ccw clipping edge, a
    // Positive value means that the point lies on the interior side.
    const auto cross = [](const Vec2& a, const Vec2& b) {
        return a(0) * b(1) - a(1) * b(0);
    };

    // Clip the current polygon successively against every edge of the clipping
    // Polygon
    for (std::size_t edge_id = 0; edge_id < clip.size(); ++edge_id) {
        // Once no points remain, later clipping edges cannot create a polygon
        if (polygon_size == 0) {
            return {};
        }

        // Define the current clipping edge by its start point and direction
        const Vec2 edge_start = clip[edge_id];
        const Vec2 edge       = clip[(edge_id + 1) % clip.size()] - edge_start;

        std::size_t clipped_size = 0;

        // Begin with the closing edge from the final current polygon point to
        // Its first point
        Vec2      previous        = polygon[polygon_size - 1];
        Precision previous_side   = cross(edge, previous - edge_start);
        bool      previous_inside = previous_side >= -tolerance;

        // Inspect every edge of the current polygon and emit the appropriate
        // Intersection and endpoint according to Sutherland-Hodgman
        for (std::size_t point_id = 0; point_id < polygon_size; ++point_id) {
            const Vec2      current        = polygon[point_id];
            const Precision current_side   = cross(edge, current - edge_start);
            const bool      current_inside = current_side >= -tolerance;

            // When the subject edge crosses the clipping boundary, insert the
            // Intersection before handling the current endpoint
            if (previous_inside != current_inside) {
                const Precision denominator = previous_side - current_side;

                // A sufficiently small denominator indicates an approximately
                // Parallel configuration for which no stable intersection
                // Should be generated
                if (std::abs(denominator) > tolerance) {
                    // Interpolate along the subject edge from previous to
                    // Current. the tolerance shifts the clipping boundary
                    // Slightly outward to retain boundary points reliably.
                    const Precision alpha = (previous_side + tolerance) / denominator;

                    logging::error(clipped_size < capacity,
                        "SurfacePolygon intersection exceeded its fixed result capacity");

                    clipped[clipped_size++] = previous + alpha * (current - previous);
                }
            }

            // Inside endpoints remain part of the clipped polygon, while
            // Outside endpoints are discarded
            if (current_inside) {
                logging::error(clipped_size < capacity,
                    "SurfacePolygon intersection exceeded its fixed result capacity");

                clipped[clipped_size++] = current;
            }

            // Continue with the next subject edge
            previous        = current;
            previous_side   = current_side;
            previous_inside = current_inside;
        }

        // Use the generated polygon as input for the next clipping edge and
        // Switch the output to the other temporary buffer
        polygon      = clipped;
        polygon_size = clipped_size;
        clipped      = clipped == buffer_a.data() ? buffer_b.data() : buffer_a.data();
    }

    // Transfer the final temporary polygon into the fixed-capacity result type
    SurfacePolygon<capacity> result;

    for (std::size_t i = 0; i < polygon_size; ++i) {
        result.push_back(polygon[i]);
    }

    return result;
}

} // namespace fem
