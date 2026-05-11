"""Small 3D vector helpers."""

from __future__ import annotations

from math import sqrt
from typing import Iterable

Point3 = tuple[float, float, float]


def point3(value: Iterable[float]) -> Point3:
    values = tuple(float(item) for item in value)
    if len(values) == 2:
        return (values[0], values[1], 0.0)
    if len(values) != 3:
        raise ValueError("point must have two or three coordinates")
    return (values[0], values[1], values[2])


def add(a: Point3, b: Point3) -> Point3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Point3, b: Point3) -> Point3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(a: Point3, value: float) -> Point3:
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a: Point3, b: Point3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Point3, b: Point3) -> Point3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(a: Point3) -> float:
    return sqrt(dot(a, a))


def distance(a: Point3, b: Point3) -> float:
    return norm(sub(a, b))


def normalize(a: Point3) -> Point3:
    length = norm(a)
    if length <= 0.0:
        raise ValueError("cannot normalize zero vector")
    return scale(a, 1.0 / length)


def plane_basis(normal: Iterable[float], reference: Iterable[float] | None = None) -> tuple[Point3, Point3, Point3]:
    n = normalize(point3(normal))
    ref = point3(reference) if reference is not None else _default_reference(n)
    ref_in_plane = sub(ref, scale(n, dot(ref, n)))
    if norm(ref_in_plane) <= 1.0e-12:
        ref_in_plane = _default_reference(n)
        ref_in_plane = sub(ref_in_plane, scale(n, dot(ref_in_plane, n)))
    u = normalize(ref_in_plane)
    v = normalize(cross(n, u))
    return u, v, n


def project_to_plane(point: Point3, origin: Point3, u: Point3, v: Point3) -> tuple[float, float]:
    rel = sub(point, origin)
    return (dot(rel, u), dot(rel, v))


def lerp(a: Point3, b: Point3, t: float) -> Point3:
    return add(scale(a, 1.0 - t), scale(b, t))


def centroid(points: Iterable[Point3]) -> Point3:
    items = tuple(points)
    if not items:
        raise ValueError("cannot compute centroid of no points")
    inv = 1.0 / len(items)
    return (
        sum(point[0] for point in items) * inv,
        sum(point[1] for point in items) * inv,
        sum(point[2] for point in items) * inv,
    )


def _default_reference(normal: Point3) -> Point3:
    x_axis = (1.0, 0.0, 0.0)
    y_axis = (0.0, 1.0, 0.0)
    return x_axis if abs(dot(normal, x_axis)) < 0.9 else y_axis
