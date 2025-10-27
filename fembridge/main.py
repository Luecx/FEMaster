# main.py
from __future__ import annotations

from pathlib import Path

if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    __package__ = "fembridge"

# core
from .core import Model, Runner

# nodes / elements
from .nodes.node import Node
from .elements.b33 import B33

# sets
from .sets.elementset import ElementSet

# materials (class-based composition)
from .materials import (
    Material,
    ElasticityIsotropic,
    Density,
)

# beam profiles
from .profiles.rect_profile import RectProfile

# sections
from .sections.beam_section import BeamSection

# supports & loads
from .supports.support import Support
from .supports.support_collector import SupportCollector
from .loads import CLoad, LoadCollector
from .steps import LinearStaticStep, EigenfrequencyStep, LinearBucklingStep


def build_beam_model() -> tuple[Model, SupportCollector, LoadCollector]:
    m = Model("Beam3Elements")

    # --- NODES: 4 nodes along X (0 -> 3) ---
    n0 = m.add_node(Node(None, 0.0, 0.0, 0.0))
    n1 = m.add_node(Node(None, 1.0, 0.0, 0.0))
    n2 = m.add_node(Node(None, 2.0, 0.0, 0.0))
    n3 = m.add_node(Node(None, 3.0, 0.0, 0.0))

    # --- ELEMENTS: 3 beam elements (B33) ---
    e0 = m.add_element(B33(None, [n0, n1]))
    e1 = m.add_element(B33(None, [n1, n2]))
    e2 = m.add_element(B33(None, [n2, n3]))

    # --- ELEMENT SET for section assignment ---
    beam_elset = m.add_elementset(ElementSet("BEAMS", [e0, e1, e2]))

    # --- MATERIAL ---
    steel = m.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e9, 0.30))
        .set_density(Density(7850.0))
    )

    # --- PROFILE ---
    rect = RectProfile("R40x20", b=0.040, h=0.020)
    m.add_profile(rect)  # ensure it gets emitted once

    # --- SECTION ---
    m.add_section(BeamSection(profile=rect, material=steel, elset=beam_elset, n1=(0.0, 1.0, 0.0)))

    # --- SUPPORTS ---
    # left node fixed in all 6 DOFs
    sup_left = Support(n0, (0, 0, 0, 0, 0, 0))

    sc = SupportCollector("SUPPORTS_LEFT")
    sc.add(sup_left)
    m.add_supportcollector(sc)

    # --- LOADS ---
    # concentrated load at the right node: downward 1 kN in Y
    lc = LoadCollector("LOAD_RIGHT", [CLoad(n3, (0.0, -1000.0, 0.0))])
    m.add_loadcollector(lc)

    return m, sc, lc


if __name__ == "__main__":
    model, support_collector, load_collector = build_beam_model()
    model.add_step(
        LinearStaticStep(
            load_collectors=[load_collector],
            support_collectors=[support_collector],
        )
    )
    model.add_step(
        EigenfrequencyStep(
            num_eigenvalues=5,
            support_collectors=[support_collector],
        )
    )

    print(model.to_femaster())
    print()
    print(model.steps.to_femaster())

    engine_path = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
    # runner.set_option(Runner.Option.NO_TEMP_FILES)
    solution = runner.run()
    print()
    print(solution)

    import numpy as np
    print(np.sqrt(np.asarray(solution.steps[1].fields["EIGENVALUES"])) / (2 * np.pi))

    for frame in solution.steps[1].frames:
        print(f"Mode {frame.index + 1} frequencies:")
        print(frame.fields["MODE_SHAPE"])