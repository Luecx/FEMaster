"""Smoke tests for the object-based public API."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from femaster_api import FEMasterWriter, Model, ResultReader
from femaster_api.model import (
    Element,
    ElementSet,
    IsotropicElasticity,
    LinearStaticLoadcase,
    LoadCollector,
    Material,
    NodalForce,
    Node,
    NodeSet,
    ShellSection,
    Support,
    SupportCollector,
)
from femaster_api.model.elements import S4
from femaster_api.tools.tovtk import build_vtk_mesh, write_vtk


class SmokeTest(unittest.TestCase):
    def test_deck_writer_and_vtu_export(self) -> None:
        model = _cantilever_model()

        diagnostics = model.validate()
        self.assertTrue(diagnostics.ok, diagnostics.errors)

        deck = FEMasterWriter(model).render()
        self.assertIn("*NODE", deck)
        self.assertIn("1, 0, 0, 0", deck)
        self.assertIn("*ELEMENT, TYPE=S4", deck)
        self.assertIn("*LOADCASE, TYPE=LINEARSTATIC, NAME=LC1", deck)

        mesh = build_vtk_mesh(model)
        self.assertEqual(len(mesh.points), 4)
        self.assertEqual(mesh.connectivity, (0, 1, 2, 3))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cantilever.vtk"
            write_vtk(model, path)
            self.assertIn("DATASET UNSTRUCTURED_GRID", path.read_text(encoding="utf-8"))

    def test_result_reader_frame_structure(self) -> None:
        result = ResultReader().parse(
            """
            LOADCASE 1 LC1
            FRAME TIME=0.25
            FIELD, NAME=DISPLACEMENT, ROWS=1, COLS=4
            1, 0.0, 0.0, -1.0
            END FIELD
            """
        )

        frame = result.loadcases["LC1"].frames[0]
        self.assertEqual(frame.time, 0.25)
        self.assertEqual(frame.fields["DISPLACEMENT"].row(0), (1.0, 0.0, 0.0, -1.0))


def _cantilever_model() -> Model:
    model = Model("cantilever")

    n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
    n2 = model.nodes.add(Node(1.0, 0.0, 0.0))
    n3 = model.nodes.add(Node(1.0, 1.0, 0.0))
    n4 = model.nodes.add(Node(0.0, 1.0, 0.0))

    plate = model.elements.add(Element(S4, [n1, n2, n3, n4]))

    fixed = model.node_sets.add(NodeSet("FIXED", [n1, n4]))
    loaded = model.node_sets.add(NodeSet("LOADED", [n2, n3]))
    plate_set = model.element_sets.add(ElementSet("PLATE_ELEMENTS", [plate]))

    steel = model.materials.add(
        Material(
            "STEEL",
            elasticity=IsotropicElasticity(
                youngs_modulus=210000.0,
                poisson_ratio=0.3,
            ),
        )
    )

    model.sections.add(
        ShellSection(
            "PLATE",
            material=steel,
            element_set=plate_set,
            thickness=1.0,
        )
    )

    bcs = model.support_collectors.add(
        SupportCollector(
            "BCS",
            [
                Support(fixed, ux=0.0, uy=0.0, uz=0.0),
            ],
        )
    )

    loads = model.load_collectors.add(
        LoadCollector(
            "LOADS",
            [
                NodalForce(loaded, fz=-1000.0),
            ],
        )
    )

    model.loadcases.add(
        LinearStaticLoadcase(
            "LC1",
            supports=[bcs],
            loads=[loads],
        )
    )

    return model


if __name__ == "__main__":
    unittest.main()
