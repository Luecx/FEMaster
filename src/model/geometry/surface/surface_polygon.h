/**
 * @file surface_polygon.h
 * @brief Declares a fixed-capacity polygon for surface operations in 2D.
 *
 * `SurfacePolygon` stores a small number of two-dimensional points without
 * dynamic memory allocation. It provides signed and unsigned area evaluation,
 * a vertex-average point, orientation handling and convex polygon intersection
 * using the Sutherland-Hodgman clipping algorithm.
 *
 * Polygon vertices used for intersection must be ordered counter-clockwise.
 * The clipping polygon and subject polygon are assumed to be convex.
 * The template capacity is fixed at compile time; intersection results use a
 * capacity equal to the sum of the two input capacities.
 *
 * @see surface_polygon.inl
 */

#pragma once

#include "../../../core/types_eig.h"

#include <array>
#include <cstddef>
#include <initializer_list>

namespace fem {

/**
 * @class SurfacePolygon
 * @brief Fixed-capacity two-dimensional polygon without dynamic allocation.
 *
 * @tparam MP Maximum number of vertices stored by the polygon.
 */
template<std::size_t MP>
class SurfacePolygon {
    // Fixed-capacity vertex storage and number of active vertices
    std::array<Vec2, MP> points_{};
    std::size_t          size_ = 0;
public:
    // Construction
    SurfacePolygon() = default;
    SurfacePolygon(std::initializer_list<Vec2> points);

    // Number of active vertices currently stored in the fixed-capacity array.
    // Unused array entries do not belong to the polygon and are ignored by all
    // geometric operations.
    [[nodiscard]] std::size_t size() const;

    // Checked access to vertices in their boundary order. The mutable overload
    // allows callers to update coordinates without changing the polygon size.
    [[nodiscard]] const Vec2& operator[](std::size_t index) const;
    [[nodiscard]] Vec2&       operator[](std::size_t index);

    // Append a vertex while preserving the existing boundary order. The
    // compile-time capacity is enforced by the implementation.
    void push_back(const Vec2& point);

    // Reverse the boundary orientation while retaining the first vertex as a
    // deterministic anchor for subsequent clipping and triangulation.
    void flip();

    // Signed and unsigned polygon area together with the arithmetic mean of
    // the active vertices. `area2()` also provides the orientation sign.
    [[nodiscard]] Precision area2() const;
    [[nodiscard]] Precision area() const;

    [[nodiscard]] Vec2 mid() const;

    // Determine whether the signed area indicates counter-clockwise vertex
    // ordering, as required by the clipping algorithm.
    [[nodiscard]] bool is_ccw() const;

    // Intersect two convex polygons using fixed-capacity Sutherland-Hodgman
    // clipping. The result capacity is the sum of both input capacities.
    template<std::size_t MPO>
    [[nodiscard]] SurfacePolygon<MP + MPO> intersection(const SurfacePolygon<MPO>& clip_polygon) const;

};

} // namespace fem

#include "surface_polygon.inl"
