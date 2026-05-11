"""Geometry helper entities for model generation and mesh tagging."""

from .base import Curve, GeometryEntity, PointEntity, SurfaceEntity, VolumeEntity
from .curves import Arc, Circle, JoinedCurve, Line, Polyline, Spline
from .meshing import GmshMesher, mesh_with_gmsh
from .surfaces import Face
from .tagging import apply_geometry_tags
from .vector import Point3, distance
from .volumes import Volume

__all__ = [
    "Arc",
    "Circle",
    "Curve",
    "Face",
    "GeometryEntity",
    "GmshMesher",
    "JoinedCurve",
    "Line",
    "Point3",
    "PointEntity",
    "Polyline",
    "Spline",
    "SurfaceEntity",
    "VolumeEntity",
    "Volume",
    "apply_geometry_tags",
    "distance",
    "mesh_with_gmsh",
]
