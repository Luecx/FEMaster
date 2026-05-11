"""Build Gmsh geometry from geometry entities."""

from __future__ import annotations

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.curves import Arc, Circle, JoinedCurve, Line, Polyline, Spline
from femaster_api.tools.geometry.surfaces import Face
from femaster_api.tools.geometry.vector import Point3
from femaster_api.tools.geometry.volumes import Volume


class GmshModelBuilder:
    def __init__(self, gmsh, mesh_size: float) -> None:
        self.gmsh = gmsh
        self.mesh_size = mesh_size
        self._points: dict[tuple[int, int, int], int] = {}
        self._lines: dict[tuple[tuple[int, int, int], tuple[int, int, int]], tuple[int, int]] = {}

    def add_face(self, face: Face, *, recombine_transfinite: bool) -> int:
        outer = self._add_loop(face.boundary)
        holes = [self._add_loop(hole) for hole in face.holes]
        surface = self.gmsh.model.geo.addPlaneSurface([outer, *holes])
        if face.is_transfinite:
            for curve in face.boundary:
                for tag in self.add_curve(curve):
                    self.gmsh.model.geo.mesh.setTransfiniteCurve(abs(tag), curve.subdivisions + 1)
            self.gmsh.model.geo.mesh.setTransfiniteSurface(surface)
            if recombine_transfinite:
                self.gmsh.model.geo.mesh.setRecombine(2, surface)
        return surface

    def add_volume(self, volume: Volume, *, recombine_transfinite: bool) -> int:
        surfaces = [self.add_face(face, recombine_transfinite=recombine_transfinite) for face in volume.faces]
        shell = self.gmsh.model.geo.addSurfaceLoop(surfaces)
        return self.gmsh.model.geo.addVolume([shell])

    def add_curve(self, curve: Curve) -> tuple[int, ...]:
        if isinstance(curve, Line):
            return (self._line(curve.a, curve.b, curve.subdivisions),)
        if isinstance(curve, Polyline):
            tags = []
            for a, b in zip(curve.points, curve.points[1:]):
                tags.append(self.gmsh.model.geo.addLine(self._point(a), self._point(b)))
            _seed_curve_tags(self.gmsh, tags, curve.subdivisions)
            return tuple(tags)
        if isinstance(curve, JoinedCurve):
            return tuple(tag for item in curve.curves for tag in self.add_curve(item))
        if isinstance(curve, Arc):
            tag = self.gmsh.model.geo.addCircleArc(self._point(curve.start), self._point(curve.center), self._point(curve.end))
            self.gmsh.model.geo.mesh.setTransfiniteCurve(tag, curve.subdivisions + 1)
            return (tag,)
        if isinstance(curve, Circle):
            return self._add_circle(curve)
        if isinstance(curve, Spline):
            points = [self._point(point) for point in curve.sample(curve.subdivisions + 1)]
            tag = self.gmsh.model.geo.addSpline(points)
            self.gmsh.model.geo.mesh.setTransfiniteCurve(tag, curve.subdivisions + 1)
            return (tag,)
        raise TypeError(f"unsupported curve: {type(curve).__name__}")

    def _add_loop(self, curves: tuple[Curve, ...]) -> int:
        tags = tuple(tag for curve in curves for tag in self.add_curve(curve))
        return self.gmsh.model.geo.addCurveLoop(list(tags))

    def _add_circle(self, circle: Circle) -> tuple[int, ...]:
        center = self._point(circle.center)
        samples = [circle.point_at(i / 4.0) for i in range(4)]
        points = [self._point(point) for point in samples]
        tags = []
        for start, end in zip(points, (*points[1:], points[0])):
            tags.append(self.gmsh.model.geo.addCircleArc(start, center, end))
        _seed_curve_tags(self.gmsh, tags, circle.subdivisions)
        return tuple(tags)

    def _point(self, point: Point3) -> int:
        key = _point_key(point)
        existing = self._points.get(key)
        if existing is not None:
            return existing
        tag = self.gmsh.model.geo.addPoint(point[0], point[1], point[2], self.mesh_size)
        self._points[key] = tag
        return tag

    def _line(self, start: Point3, end: Point3, subdivisions: int) -> int:
        start_key = _point_key(start)
        end_key = _point_key(end)
        key = (start_key, end_key)
        reverse_key = (end_key, start_key)
        if key in self._lines:
            tag, existing_subdivisions = self._lines[key]
            _validate_same_line_seed(existing_subdivisions, subdivisions, start_key, end_key)
            return tag
        if reverse_key in self._lines:
            tag, existing_subdivisions = self._lines[reverse_key]
            _validate_same_line_seed(existing_subdivisions, subdivisions, start_key, end_key)
            return -tag
        tag = self.gmsh.model.geo.addLine(self._point(start), self._point(end))
        self.gmsh.model.geo.mesh.setTransfiniteCurve(tag, subdivisions + 1)
        self._lines[key] = (tag, subdivisions)
        return tag


def _seed_curve_tags(gmsh, tags: list[int], subdivisions: int) -> None:
    per_curve = max(1, subdivisions // max(1, len(tags)))
    for tag in tags:
        gmsh.model.geo.mesh.setTransfiniteCurve(abs(tag), per_curve + 1)


def _point_key(point: Point3) -> tuple[int, int, int]:
    return (round(point[0] * 1.0e9), round(point[1] * 1.0e9), round(point[2] * 1.0e9))


def _validate_same_line_seed(
    existing_subdivisions: int,
    subdivisions: int,
    start_key: tuple[int, int, int],
    end_key: tuple[int, int, int],
) -> None:
    if existing_subdivisions != subdivisions:
        raise ValueError(
            "duplicate geometry line has conflicting seeds: "
            f"{existing_subdivisions} != {subdivisions} for edge {_format_point_key(start_key)} -> {_format_point_key(end_key)}"
        )


def _format_point_key(point: tuple[int, int, int]) -> tuple[float, float, float]:
    return (point[0] / 1.0e9, point[1] / 1.0e9, point[2] / 1.0e9)
