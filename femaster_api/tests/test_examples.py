"""Regression tests for the simplified solid beam example."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from femaster_api import ResultReader
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
    SolidSection,
    Support,
    SupportCollector,
)
from femaster_api.model.elements import C3D8
from femaster_api.model.model import Model
from femaster_api.tools.tovtk import build_vtk_mesh, write_vtk


class SolidBeamExampleTest(unittest.TestCase):
    def test_example_script_runs_with_fake_femaster(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(__file__).resolve().parents[1]
            fake = Path(directory) / "fake_femaster.py"
            fake.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env python3
                    import sys
                    from pathlib import Path

                    output = Path(sys.argv[sys.argv.index("--output") + 1])
                    rows = "\\n".join("0.0, 0.0, -1.0" for _ in range(29))
                    output.write_text(
                        f"LOADCASE 1 LC_BENDING\\n"
                        f"FIELD, NAME=DISPLACEMENT, ROWS=29, COLS=3\\n"
                        f"{rows}\\n"
                        f"END FIELD\\n",
                        encoding="utf-8",
                    )
                    print("fake femaster completed")
                    """
                ),
                encoding="utf-8",
            )
            fake.chmod(0o755)

            env = os.environ.copy()
            env["FEMASTER_EXECUTABLE"] = str(fake)
            env["PYTHONPATH"] = f"{root}{os.pathsep}{env.get('PYTHONPATH', '')}"

            completed = subprocess.run(
                [sys.executable, str(root / "examples" / "solid_beam.py")],
                cwd=directory,
                env=env,
                text=True,
                capture_output=True,
                check=False,
            )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Result summary:", completed.stdout)
        self.assertIn("DISPLACEMENT", completed.stdout)
        self.assertIn("Wrote VTK:", completed.stdout)

    def test_vtk_writer_includes_point_results(self) -> None:
        model = _solid_beam_model()
        rows = "\n".join("0.0, 0.0, -1.0" for _ in range(len(model.nodes) + 1))
        result = ResultReader().parse(
            f"""
            LOADCASE 1 LC_BENDING
            FIELD, NAME=DISPLACEMENT, ROWS={len(model.nodes) + 1}, COLS=3
            {rows}
            END FIELD
            """
        )

        with tempfile.TemporaryDirectory() as directory:
            vtk_path = Path(directory) / "solid_beam.vtk"
            write_vtk(model, vtk_path, result=result)
            text = vtk_path.read_text(encoding="utf-8")

        self.assertIn("POINT_DATA 28", text)
        self.assertIn("DISPLACEMENT", text)

    def test_mesh_builder_handles_solid_beam(self) -> None:
        mesh = build_vtk_mesh(_solid_beam_model())

        self.assertEqual(len(mesh.points), 28)
        self.assertEqual(len(mesh.cell_types), 6)
        self.assertEqual(set(mesh.cell_types), {12})


def _solid_beam_model() -> Model:
    model = Model("solid_beam")

    nodes: dict[tuple[int, int, int], Node] = {}
    for ix in range(7):
        for iy in range(2):
            for iz in range(2):
                nodes[(ix, iy, iz)] = model.nodes.add(Node(float(ix), float(iy), float(iz)))

    elements = []
    for ix in range(6):
        elements.append(
            model.elements.add(
                Element(
                    C3D8,
                    [
                        nodes[(ix, 0, 0)],
                        nodes[(ix + 1, 0, 0)],
                        nodes[(ix + 1, 1, 0)],
                        nodes[(ix, 1, 0)],
                        nodes[(ix, 0, 1)],
                        nodes[(ix + 1, 0, 1)],
                        nodes[(ix + 1, 1, 1)],
                        nodes[(ix, 1, 1)],
                    ],
                )
            )
        )

    fixed = model.node_sets.add(NodeSet("FIXED", [node for (ix, _, _), node in nodes.items() if ix == 0]))
    loaded = model.node_sets.add(NodeSet("LOADED_EDGE", [node for (ix, _, iz), node in nodes.items() if ix == 6 and iz == 1]))
    solids = model.element_sets.add(ElementSet("BEAM_SOLIDS", elements))

    steel = model.materials.add(
        Material(
            "STEEL",
            elasticity=IsotropicElasticity(youngs_modulus=210000.0, poisson_ratio=0.3),
        )
    )
    model.sections.add(SolidSection("BEAM_SECTION", material=steel, element_set=solids))

    supports = model.support_collectors.add(
        SupportCollector("SUPPORTS", [Support(fixed, ux=0.0, uy=0.0, uz=0.0)])
    )
    loads = model.load_collectors.add(LoadCollector("LOADS", [NodalForce(loaded, fz=-500.0)]))
    model.loadcases.add(LinearStaticLoadcase("LC_BENDING", supports=[supports], loads=[loads]))
    return model


if __name__ == "__main__":
    unittest.main()
