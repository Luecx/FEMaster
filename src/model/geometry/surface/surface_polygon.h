/**
 * @file surface_polygon.h
 * @brief Declares a fixed-capacity polygon for surface operations in 2D.
 *
 * `SurfacePolygon` stores a small number of two-dimensional points without
 * dynamic memory allocation. It provides polygon area evaluation and convex
 * polygon intersection using the Sutherland-Hodgman clipping algorithm.
 *
 * Polygon vertices used for intersection must be ordered counter-clockwise.
 * The clipping polygon and subject polygon are assumed to be convex.
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
    // points and current amount of points actually stored/used
    std::array<Vec2, MP> points_{};
    std::size_t          size_ = 0;
public:
    // constructions
    SurfacePolygon() = default;
    SurfacePolygon(std::initializer_list<Vec2> points);

    // basic functions for size
    [[nodiscard]] std::size_t size() const;

    // accessor functions
    [[nodiscard]] const Vec2& operator[](std::size_t index) const;
    [[nodiscard]] Vec2&       operator[](std::size_t index);
    // add a new point
    void push_back(const Vec2& point);

    // reverse the polygon orientation while retaining the first point
    void flip();

    // area functions, area2 is just twice the area and has a sign convention and may be negative if not ccw
    [[nodiscard]] Precision area2() const;
    [[nodiscard]] Precision area() const;

    // mid points of the polygon
    [[nodiscard]] Vec2      mid() const;

    // check if its counter clock wise, if
    [[nodiscard]] bool is_ccw() const;

    // intersection method
    template<std::size_t MPO>
    [[nodiscard]] SurfacePolygon<MP + MPO> intersection(const SurfacePolygon<MPO>& clip_polygon) const;

};

} // namespace fem

#include "surface_polygon.inl"