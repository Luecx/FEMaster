/**
 * @file surface_polygon.inl
 * @brief Implements the fixed-capacity polygon for surface operations in 2D.
 *
 * @see surface_polygon.h
 */

#pragma once

#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"

#include <cmath>
#include <utility>

namespace fem {

template<std::size_t MP>
SurfacePolygon<MP>::SurfacePolygon(std::initializer_list<Vec2> points) {
    // the number of supplied points must fit into the fixed-capacity storage
    logging::error(points.size() <= MP,
        "SurfacePolygon initializer contains more points than its fixed capacity");

    // copy the supplied points into the internal array and update the actual polygon size
    for (const Vec2& point : points) {
        points_[size_++] = point;
    }
}

template<std::size_t MP>
std::size_t SurfacePolygon<MP>::size() const {
    return size_;
}

template<std::size_t MP>
const Vec2& SurfacePolygon<MP>::operator[](std::size_t index) const {
    // only the first size_ entries of points_ belong to the actual polygon
    logging::error(index < size_,
        "SurfacePolygon point index is out of bounds");

    return points_[index];
}

template<std::size_t MP>
Vec2& SurfacePolygon<MP>::operator[](std::size_t index) {
    // only the first size_ entries of points_ belong to the actual polygon
    logging::error(index < size_,
        "SurfacePolygon point index is out of bounds");

    return points_[index];
}

template<std::size_t MP>
void SurfacePolygon<MP>::push_back(const Vec2& point) {
    // SurfacePolygon has fixed storage and must never allocate dynamically
    logging::error(size_ < MP,
        "SurfacePolygon fixed capacity exceeded while adding a point");

    points_[size_++] = point;
}

template<std::size_t MP>
void SurfacePolygon<MP>::flip() {
    // polygons containing fewer than three points do not have a meaningful orientation
    if (size_ < 3) {
        return;
    }

    // reverse the winding while retaining the first point. keeping the first
    // point fixed makes the operation deterministic with respect to the polygon
    // anchor while still changing between clockwise and counter-clockwise order.
    for (std::size_t left = 1, right = size_ - 1; left < right; ++left, --right) {
        std::swap(points_[left], points_[right]);
    }
}

template<std::size_t MP>
Precision SurfacePolygon<MP>::area2() const {
    Precision result = Precision(0);

    // evaluate the shoelace formula over all polygon edges. the modulo closes
    // the polygon by connecting the final point back to the first point.
    for (std::size_t i = 0; i < size_; ++i) {
        const Vec2& a = points_[i];
        const Vec2& b = points_[(i + 1) % size_];

        result += a(0) * b(1) - a(1) * b(0);
    }

    // the sign contains the polygon orientation:
    // positive for counter-clockwise and negative for clockwise ordering
    return result;
}

template<std::size_t MP>
Vec2 SurfacePolygon<MP>::mid() const {
    // an empty polygon does not provide a meaningful average point
    logging::error(size_ > 0,
        "SurfacePolygon cannot compute the midpoint of an empty polygon");
    
    Vec2 result{0,0};

    // just sum over all points and then divide by amount of points
    for (std::size_t i = 0; i < size_; ++i) {
        result += points_[i];
    }

    // divide to adjust for summation above
    return result / size_;
}

template<std::size_t MP>
Precision SurfacePolygon<MP>::area() const {
    // area2 is signed, while the physical polygon area is always non-negative
    return Precision(0.5) * std::abs(area2());
}

template<std::size_t MP>
bool SurfacePolygon<MP>::is_ccw() const {
    // a non-negative signed area indicates counter-clockwise vertex ordering
    return area2() >= Precision(0);
}

template<std::size_t MP>
template<std::size_t MPO>
SurfacePolygon<MP + MPO> SurfacePolygon<MP>::intersection(
    const SurfacePolygon<MPO>& clip_polygon
) const {
    // the intersection of two convex polygons with MP and MPO vertices can
    // contain at most MP + MPO vertices
    constexpr std::size_t capacity  = MP + MPO;
    constexpr Precision   tolerance = Precision(1e-12);

    // a valid polygonal intersection requires both inputs to contain at least
    // three points
    if (size_ < 3 || clip_polygon.size() < 3) {
        return {};
    }

    // Sutherland-Hodgman interprets the left side of every clipping edge as
    // the interior half-plane. therefore the clipping polygon has to be ccw.
    // we normalize both polygons so the resulting intersection is also ccw.
    SurfacePolygon<MP>  subject = *this;
    SurfacePolygon<MPO> clip    = clip_polygon;

    if (!subject.is_ccw()) {
        subject.flip();
    }

    if (!clip.is_ccw()) {
        clip.flip();
    }

    // clipping is performed repeatedly. two fixed-capacity buffers allow us
    // to alternate between input and output without dynamic allocation.
    std::array<Vec2, capacity> buffer_a{};
    std::array<Vec2, capacity> buffer_b{};

    // the first clipping pass reads directly from the normalized subject
    // polygon. later passes read from the temporary clipping buffers.
    const Vec2* polygon      = subject.points_.data();
    Vec2*       clipped      = buffer_a.data();
    std::size_t polygon_size = subject.size_;

    // signed two-dimensional cross product. for a ccw clipping edge, a
    // positive value means that the point lies on the interior side.
    const auto cross = [](const Vec2& a, const Vec2& b) {
        return a(0) * b(1) - a(1) * b(0);
    };

    // clip the current polygon successively against every edge of the clipping
    // polygon
    for (std::size_t edge_id = 0; edge_id < clip.size(); ++edge_id) {
        // once no points remain, later clipping edges cannot create a polygon
        if (polygon_size == 0) {
            return {};
        }

        // define the current clipping edge by its start point and direction
        const Vec2 edge_start = clip[edge_id];
        const Vec2 edge       = clip[(edge_id + 1) % clip.size()] - edge_start;

        std::size_t clipped_size = 0;

        // begin with the closing edge from the final current polygon point to
        // its first point
        Vec2      previous        = polygon[polygon_size - 1];
        Precision previous_side   = cross(edge, previous - edge_start);
        bool      previous_inside = previous_side >= -tolerance;

        // inspect every edge of the current polygon and emit the appropriate
        // intersection and endpoint according to Sutherland-Hodgman
        for (std::size_t point_id = 0; point_id < polygon_size; ++point_id) {
            const Vec2      current        = polygon[point_id];
            const Precision current_side   = cross(edge, current - edge_start);
            const bool      current_inside = current_side >= -tolerance;

            // when the subject edge crosses the clipping boundary, insert the
            // intersection before handling the current endpoint
            if (previous_inside != current_inside) {
                const Precision denominator = previous_side - current_side;

                // a sufficiently small denominator indicates an approximately
                // parallel configuration for which no stable intersection
                // should be generated
                if (std::abs(denominator) > tolerance) {
                    // interpolate along the subject edge from previous to
                    // current. the tolerance shifts the clipping boundary
                    // slightly outward to retain boundary points reliably.
                    const Precision alpha = (previous_side + tolerance) / denominator;

                    logging::error(clipped_size < capacity,
                        "SurfacePolygon intersection exceeded its fixed result capacity");

                    clipped[clipped_size++] = previous + alpha * (current - previous);
                }
            }

            // inside endpoints remain part of the clipped polygon, while
            // outside endpoints are discarded
            if (current_inside) {
                logging::error(clipped_size < capacity,
                    "SurfacePolygon intersection exceeded its fixed result capacity");

                clipped[clipped_size++] = current;
            }

            // continue with the next subject edge
            previous        = current;
            previous_side   = current_side;
            previous_inside = current_inside;
        }

        // use the generated polygon as input for the next clipping edge and
        // switch the output to the other temporary buffer
        polygon      = clipped;
        polygon_size = clipped_size;
        clipped      = clipped == buffer_a.data() ? buffer_b.data() : buffer_a.data();
    }

    // transfer the final temporary polygon into the fixed-capacity result type
    SurfacePolygon<capacity> result;

    for (std::size_t i = 0; i < polygon_size; ++i) {
        result.push_back(polygon[i]);
    }

    return result;
}

} // namespace fem