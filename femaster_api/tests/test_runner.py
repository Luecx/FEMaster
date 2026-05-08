"""Tests for high-level FEMaster execution."""

from __future__ import annotations

import os
import tempfile
import textwrap
import unittest
from pathlib import Path

from femaster_api import FEMaster, Result
from femaster_api.backend import find_femaster_executable
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


def build_model() -> Model:
    model = Model("runner_test")

    n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
    n2 = model.nodes.add(Node(1.0, 0.0, 0.0))
    n3 = model.nodes.add(Node(1.0, 1.0, 0.0))
    n4 = model.nodes.add(Node(0.0, 1.0, 0.0))
    n5 = model.nodes.add(Node(0.0, 0.0, 1.0))
    n6 = model.nodes.add(Node(1.0, 0.0, 1.0))
    n7 = model.nodes.add(Node(1.0, 1.0, 1.0))
    n8 = model.nodes.add(Node(0.0, 1.0, 1.0))

    element = model.elements.add(Element(C3D8, [n1, n2, n3, n4, n5, n6, n7, n8]))
    fixed = model.node_sets.add(NodeSet("FIXED", [n1, n4, n5, n8]))
    loaded = model.node_sets.add(NodeSet("LOADED", [n2, n3, n6, n7]))
    solids = model.element_sets.add(ElementSet("SOLIDS", [element]))

    steel = model.materials.add(
        Material("STEEL", elasticity=IsotropicElasticity(youngs_modulus=210000.0, poisson_ratio=0.3))
    )
    model.sections.add(SolidSection("SOLID_SECTION", material=steel, element_set=solids))

    supports = model.support_collectors.add(
        SupportCollector("SUPPORTS", [Support(fixed, ux=0.0, uy=0.0, uz=0.0)])
    )
    loads = model.load_collectors.add(LoadCollector("LOADS", [NodalForce(loaded, fz=-1.0)]))
    model.loadcases.add(LinearStaticLoadcase("LC_BENDING", supports=[supports], loads=[loads]))
    return model


class RunnerTest(unittest.TestCase):
    @unittest.skipIf(os.name == "nt", "script execution test uses a POSIX shebang")
    def test_run_returns_parsed_result(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fake = Path(directory) / "fake_femaster"
            fake.write_text(
                textwrap.dedent(
                    """\
                    #!/usr/bin/env python3
                    import sys
                    from pathlib import Path

                    output = Path(sys.argv[sys.argv.index("--output") + 1])
                    output.write_text(
                        "LOADCASE 1 LC_BENDING\\n"
                        "FIELD, NAME=DISPLACEMENT, ROWS=1, COLS=4\\n"
                        "1, 0.0, 0.0, -1.0\\n"
                        "END FIELD\\n",
                        encoding="utf-8",
                    )
                    print("fake femaster completed")
                    """
                ),
                encoding="utf-8",
            )
            fake.chmod(0o755)

            solver = FEMaster(
                build_model(),
                executable=fake,
                threads=4,
                temporary_files=True,
                redirect_stdout=True,
            )

            result = solver.run()

            self.assertIsInstance(result, Result)
            self.assertIn("LC_BENDING", result.loadcases)
            self.assertEqual(result.field("DISPLACEMENT", loadcase="LC_BENDING").row(0)[3], -1.0)
            self.assertEqual(solver.last_run.returncode, 0)
            self.assertIn("fake femaster completed", solver.last_run.stdout)

    def test_missing_executable_has_clear_error(self) -> None:
        with self.assertRaises(FileNotFoundError):
            FEMaster(build_model(), executable="/definitely/not/femaster")

    def test_executable_can_be_given_as_bin_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            bin_dir = Path(directory) / "bin"
            bin_dir.mkdir()
            fake = bin_dir / "FEMaster"
            fake.write_text("", encoding="utf-8")

            self.assertEqual(find_femaster_executable(bin_dir), fake.resolve())


if __name__ == "__main__":
    unittest.main()
