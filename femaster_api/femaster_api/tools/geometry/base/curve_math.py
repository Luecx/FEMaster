"""Curve helper functions."""

from __future__ import annotations

from math import pi

from femaster_api.tools.geometry.vector import Point3, add, distance, dot, lerp, scale, sub


def point_on_segment(p: Point3, a: Point3, b: Point3, tolerance: float) -> bool:
    ab = sub(b, a)
    ap = sub(p, a)
    length2 = dot(ab, ab)
    if length2 <= tolerance * tolerance:
        return distance(p, a) <= tolerance
    u = dot(ap, ab) / length2
    if u < -tolerance or u > 1.0 + tolerance:
        return False
    closest = add(a, scale(ab, u))
    return distance(p, closest) <= tolerance


def curve_length(curve) -> float:
    points = curve.sample(max(16, curve.subdivisions * 4 + 1))
    return sum(distance(a, b) for a, b in zip(points, points[1:]))


def segment_lengths(points: tuple[Point3, ...]) -> tuple[float, ...]:
    return tuple(distance(a, b) for a, b in zip(points, points[1:]))


def catmull_rom(points: tuple[Point3, ...], t: float) -> Point3:
    if len(points) == 2:
        return lerp(points[0], points[1], t)
    scaled_t = t * (len(points) - 1)
    i = min(int(scaled_t), len(points) - 2)
    local = scaled_t - i
    p0 = points[max(0, i - 1)]
    p1 = points[i]
    p2 = points[i + 1]
    p3 = points[min(len(points) - 1, i + 2)]
    local2 = local * local
    local3 = local2 * local
    return (
        0.5
        * (
            2 * p1[0]
            + (-p0[0] + p2[0]) * local
            + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * local2
            + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * local3
        ),
        0.5
        * (
            2 * p1[1]
            + (-p0[1] + p2[1]) * local
            + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * local2
            + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * local3
        ),
        0.5
        * (
            2 * p1[2]
            + (-p0[2] + p2[2]) * local
            + (2 * p0[2] - 5 * p1[2] + 4 * p2[2] - p3[2]) * local2
            + (-p0[2] + 3 * p1[2] - 3 * p2[2] + p3[2]) * local3
        ),
    )


def angle_between(angle: float, start: float, end: float, tolerance: float) -> bool:
    tau = 2.0 * pi
    angle = angle % tau
    start = start % tau
    end = end % tau
    if start <= end:
        return start - tolerance <= angle <= end + tolerance
    return angle >= start - tolerance or angle <= end + tolerance


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))
