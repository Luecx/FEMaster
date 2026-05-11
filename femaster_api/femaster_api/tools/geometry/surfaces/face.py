"""Planar surface geometry entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base import Curve, SurfaceEntity
from femaster_api.tools.geometry.vector import Point3, cross, distance, dot, norm, plane_basis, point3, project_to_plane, sub


@dataclass(frozen=True, slots=True, init=False)
class Face(SurfaceEntity):
    boundary: tuple[Curve, ...]
    holes: tuple[tuple[Curve, ...], ...] = ()
    tags: tuple[str, ...] = ()

    def __init__(
        self,
        boundary: Iterable[Curve],
        *,
        holes: Iterable[Iterable[Curve]] = (),
        tags: Iterable[str] = (),
    ) -> None:
        outer = tuple(boundary)
        inner = tuple(tuple(hole) for hole in holes)
        _validate_loop(outer, "boundary")
        for index, hole in enumerate(inner):
            _validate_loop(hole, f"hole {index}")
        object.__setattr__(self, "boundary", outer)
        object.__setattr__(self, "holes", inner)
        object.__setattr__(self, "tags", tuple(tags))

    @property
    def is_transfinite(self) -> bool:
        return self.transfinite_seeds is not None

    @property
    def transfinite_seeds(self) -> tuple[int, int] | None:
        return _transfinite_seeds(self.boundary, self.holes)

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        p = point3(point)
        outer = _loop_points(self.boundary)
        if any(curve.contains(p, tolerance=tolerance) for curve in self.boundary):
            return True
        origin, u, v, n = _surface_frame(outer)
        if abs(dot(sub(p, origin), n)) > tolerance:
            return False
        projected = project_to_plane(p, origin, u, v)
        if not _point_in_polygon(projected, _project_loop(outer, origin, u, v)):
            return False
        for hole in self.holes:
            if any(curve.contains(p, tolerance=tolerance) for curve in hole):
                return False
            if _point_in_polygon(projected, _project_loop(_loop_points(hole), origin, u, v)):
                return False
        return True


def _validate_loop(curves: tuple[Curve, ...], label: str) -> None:
    if len(curves) < 1:
        raise ValueError(f"{label} must contain at least one curve")
    if len(curves) == 1 and distance(curves[0].start, curves[0].end) <= 1.0e-6:
        return
    for left, right in zip(curves, (*curves[1:], curves[0])):
        if distance(left.end, right.start) > 1.0e-6:
            raise ValueError(f"{label} curves must form a closed loop")


def _transfinite_seeds(boundary: tuple[Curve, ...], holes: tuple[tuple[Curve, ...], ...], *, tolerance: float = 1.0e-6) -> tuple[int, int] | None:
    if len(boundary) != 4 or holes:
        return None
    c0, c1, c2, c3 = boundary
    if not _connected(c0, c1, tolerance):
        return None
    if not _connected(c1, c2, tolerance):
        return None
    if not _connected(c2, c3, tolerance):
        return None
    if not _connected(c3, c0, tolerance):
        return None
    if c0.subdivisions != c2.subdivisions:
        return None
    if c1.subdivisions != c3.subdivisions:
        return None
    return c0.subdivisions, c1.subdivisions


def _connected(left: Curve, right: Curve, tolerance: float) -> bool:
    return distance(left.end, right.start) <= tolerance


def _loop_points(curves: tuple[Curve, ...]) -> tuple[Point3, ...]:
    points: list[Point3] = []
    for curve in curves:
        samples = curve.sample(max(2, curve.subdivisions + 1))
        points.extend(samples[:-1])
    return tuple(points)


def _surface_frame(points: tuple[Point3, ...]) -> tuple[Point3, Point3, Point3, Point3]:
    if len(points) < 3:
        raise ValueError("surface loop requires at least three points")
    origin = points[0]
    for i in range(1, len(points) - 1):
        normal = cross(sub(points[i], origin), sub(points[i + 1], origin))
        if norm(normal) > 1.0e-12:
            u, v, n = plane_basis(normal, sub(points[i], origin))
            return origin, u, v, n
    raise ValueError("surface loop points are collinear")


def _project_loop(points: tuple[Point3, ...], origin: Point3, u: Point3, v: Point3) -> tuple[tuple[float, float], ...]:
    return tuple(project_to_plane(point, origin, u, v) for point in points)


def _point_in_polygon(point: tuple[float, float], polygon: tuple[tuple[float, float], ...]) -> bool:
    inside = False
    x, y = point
    count = len(polygon)
    if count < 3:
        return False
    j = count - 1
    for i in range(count):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        intersects = (yi > y) != (yj > y)
        if intersects:
            x_intersection = (xj - xi) * (y - yi) / (yj - yi) + xi
            if x < x_intersection:
                inside = not inside
        j = i
    return inside
