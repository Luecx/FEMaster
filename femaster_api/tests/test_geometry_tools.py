from math import pi

import pytest

from femaster_api import B33, Element, ElementTopology, Model, Node
from femaster_api.model.sets import EntityType
from femaster_api.tools.geometry import (
    Arc,
    Face,
    Circle,
    JoinedCurve,
    Line,
    PointEntity,
    Spline,
    Volume,
    apply_geometry_tags,
    mesh_with_gmsh,
)


def test_curves_containment_and_joining() -> None:
    left = Line((0.0, 0.0), (1.0, 0.0), seed=2, tags=("edge",))
    right = Arc((1.0, 1.0), 1.0, -pi / 2.0, 0.0, seed=4, tags=("arc",))
    joined = JoinedCurve((left, right), tags=("path",))

    assert left.contains((0.5, 0.0, 0.0))
    assert right.contains((2.0, 1.0, 0.0))
    assert joined.subdivisions == 6
    assert joined.contains((0.5, 0.0, 0.0))


def test_arc_and_circle_are_planar_3d_curves() -> None:
    arc = Arc((0.0, 0.0, 0.0), 1.0, 0.0, pi / 2.0, normal=(0.0, 1.0, 0.0), reference=(1.0, 0.0, 0.0))
    circle = Circle((0.0, 0.0, 0.0), 2.0, normal=(1.0, 0.0, 0.0), reference=(0.0, 1.0, 0.0))

    assert arc.point_at(1.0) == pytest.approx((0.0, 0.0, -1.0))
    assert arc.contains((0.0, 0.0, -1.0))
    assert not arc.contains((0.0, 0.1, -1.0))
    assert circle.contains((0.0, 2.0, 0.0))
    assert circle.contains((0.0, 0.0, 2.0))
    assert not circle.contains((0.1, 0.0, 2.0))


def test_spline_interpolations() -> None:
    linear = Spline(((0.0, 0.0), (1.0, 1.0)), interpolation="linear", seed=4)
    bezier = Spline(((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)), interpolation="bezier", seed=8)

    assert linear.point_at(0.5) == (0.5, 0.5, 0.0)
    assert bezier.point_at(0.0) == (0.0, 0.0, 0.0)
    assert bezier.point_at(1.0) == (1.0, 0.0, 0.0)


def test_face_holes_and_transfinite_detection() -> None:
    bottom = Line((0.0, 0.0), (1.0, 0.0), seed=4)
    right = Line((1.0, 0.0), (1.0, 1.0), seed=3)
    top = Line((1.0, 1.0), (0.0, 1.0), seed=4)
    left = Line((0.0, 1.0), (0.0, 0.0), seed=3)
    face = Face((bottom, right, top, left), tags=("plate",))

    assert face.contains((0.5, 0.5, 0.0))
    assert face.is_transfinite
    assert face.transfinite_seeds == (4, 3)


def test_face_contains_works_outside_xy_plane() -> None:
    bottom = Line((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), seed=2)
    right = Line((0.0, 1.0, 0.0), (0.0, 1.0, 1.0), seed=2)
    top = Line((0.0, 1.0, 1.0), (0.0, 0.0, 1.0), seed=2)
    left = Line((0.0, 0.0, 1.0), (0.0, 0.0, 0.0), seed=2)
    face = Face((bottom, right, top, left))

    assert face.contains((0.0, 0.5, 0.5))
    assert not face.contains((0.1, 0.5, 0.5))


def test_volume_is_constructed_from_faces() -> None:
    volume = _cube_volume(tags=("SOLID",))

    assert volume.contains((0.5, 0.5, 0.5))
    assert volume.contains((0.0, 0.5, 0.5))
    assert not volume.contains((1.5, 0.5, 0.5))


def test_geometry_tags_create_node_and_element_sets() -> None:
    model = Model("tagging")
    n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
    n2 = model.nodes.add(Node(1.0, 0.0, 0.0))
    element = model.elements.add(Element(B33, (n1, n2)))

    apply_geometry_tags(model, (PointEntity((0.0, 0.0, 0.0), tags=("LEFT",)), Line((0.0, 0.0), (1.0, 0.0), tags=("AXIS",))))

    left_nodes = model.sets.get(EntityType.NODE, "LEFT")
    axis_nodes = model.sets.get(EntityType.NODE, "AXIS")
    axis_elements = model.sets.get(EntityType.ELEMENT, "AXIS")

    assert left_nodes.members == (n1,)
    assert set(axis_nodes.members) == {n1, n2}
    assert axis_elements.members == (element,)


def test_gmsh_mesher_creates_surface_mesh_and_tag_sets() -> None:
    face = Face(
        (
            Line((0.0, 0.0), (1.0, 0.0), seed=2),
            Line((1.0, 0.0), (1.0, 1.0), seed=2),
            Line((1.0, 1.0), (0.0, 1.0), seed=2),
            Line((0.0, 1.0), (0.0, 0.0), seed=2),
        ),
        tags=("SURFACE",),
    )

    model = mesh_with_gmsh((face,), mesh_size=0.5, model_name="gmsh_test")

    assert len(model.nodes) > 0
    assert len(model.elements) > 0
    assert model.sets.get(EntityType.NODE, "SURFACE").members
    assert model.sets.get(EntityType.ELEMENT, "SURFACE").members


def test_gmsh_mesher_rejects_duplicate_line_with_conflicting_seed() -> None:
    left_face = Face(
        (
            Line((0.0, 0.0), (1.0, 0.0), seed=2),
            Line((1.0, 0.0), (1.0, 1.0), seed=2),
            Line((1.0, 1.0), (0.0, 1.0), seed=2),
            Line((0.0, 1.0), (0.0, 0.0), seed=2),
        )
    )
    right_face = Face(
        (
            Line((1.0, 0.0), (2.0, 0.0), seed=2),
            Line((2.0, 0.0), (2.0, 1.0), seed=3),
            Line((2.0, 1.0), (1.0, 1.0), seed=2),
            Line((1.0, 1.0), (1.0, 0.0), seed=3),
        )
    )

    with pytest.raises(ValueError, match="conflicting seeds"):
        mesh_with_gmsh((left_face, right_face), mesh_size=0.5, model_name="conflicting_seed_test")


def test_gmsh_mesher_creates_volume_mesh_from_faces() -> None:
    volume = _cube_volume(tags=("SOLID",))

    model = mesh_with_gmsh(
        (volume,),
        mesh_size=0.75,
        model_name="gmsh_volume_test",
        recombine_transfinite=False,
    )

    assert len(model.nodes) > 0
    assert any(element.topology is ElementTopology.TET4 for element in model.elements)
    assert model.sets.get(EntityType.NODE, "SOLID").members
    assert model.sets.get(EntityType.ELEMENT, "SOLID").members


def _cube_volume(*, tags: tuple[str, ...] = ()) -> Volume:
    p000 = (0.0, 0.0, 0.0)
    p100 = (1.0, 0.0, 0.0)
    p110 = (1.0, 1.0, 0.0)
    p010 = (0.0, 1.0, 0.0)
    p001 = (0.0, 0.0, 1.0)
    p101 = (1.0, 0.0, 1.0)
    p111 = (1.0, 1.0, 1.0)
    p011 = (0.0, 1.0, 1.0)
    faces = (
        _quad_face(p000, p100, p110, p010),
        _quad_face(p001, p011, p111, p101),
        _quad_face(p000, p001, p101, p100),
        _quad_face(p010, p110, p111, p011),
        _quad_face(p000, p010, p011, p001),
        _quad_face(p100, p101, p111, p110),
    )
    return Volume(faces, tags=tags)


def _quad_face(a, b, c, d) -> Face:
    return Face(
        (
            Line(a, b, seed=2),
            Line(b, c, seed=2),
            Line(c, d, seed=2),
            Line(d, a, seed=2),
        )
    )
