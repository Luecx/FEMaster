"""Volume entity constructed from boundary faces."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base import VolumeEntity
from femaster_api.tools.geometry.surfaces import Face
from femaster_api.tools.geometry.vector import Point3, add, centroid, cross, dot, norm, point3, scale, sub


@dataclass(frozen=True, slots=True, init=False)
class Volume(VolumeEntity):
    faces: tuple[Face, ...]
    tags: tuple[str, ...] = ()

    def __init__(self, faces: Iterable[Face], *, tags: Iterable[str] = ()) -> None:
        items = tuple(faces)
        if len(items) < 4:
            raise ValueError("volume requires at least four boundary faces")
        object.__setattr__(self, "faces", items)
        object.__setattr__(self, "tags", tuple(tags))

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        p = point3(point)
        if any(face.contains(p, tolerance=tolerance) for face in self.faces):
            return True

        intersections: list[float] = []
        ray = _ray_direction()
        for face in self.faces:
            face_center = _face_centroid(face)
            normal = _face_normal(face)
            denominator = dot(normal, ray)
            if abs(denominator) <= tolerance:
                continue
            distance_on_ray = dot(normal, sub(face_center, p)) / denominator
            if distance_on_ray <= tolerance:
                continue
            hit = add(p, scale(ray, distance_on_ray))
            if face.contains(hit, tolerance=tolerance):
                _append_unique_intersection(intersections, distance_on_ray, tolerance)
        return len(intersections) % 2 == 1


def _face_points(face: Face) -> tuple[Point3, ...]:
    points: list[Point3] = []
    for curve in face.boundary:
        points.extend(curve.sample(curve.subdivisions + 1)[:-1])
    return tuple(points)


def _face_centroid(face: Face) -> Point3:
    return centroid(_face_points(face))


def _face_normal(face: Face) -> Point3:
    points = _face_points(face)
    origin = points[0]
    for i in range(1, len(points) - 1):
        normal = cross(sub(points[i], origin), sub(points[i + 1], origin))
        length = norm(normal)
        if length > 1.0e-12:
            return (normal[0] / length, normal[1] / length, normal[2] / length)
    raise ValueError("volume face has collinear boundary points")


def _ray_direction() -> Point3:
    return (0.812286, 0.471391, 0.343724)


def _append_unique_intersection(values: list[float], value: float, tolerance: float) -> None:
    if not any(abs(existing - value) <= tolerance for existing in values):
        values.append(value)
