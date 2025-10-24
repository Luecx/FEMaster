# main.py
from __future__ import annotations

# core
from core.model import Model

# nodes / elements
from nodes.node import Node
from elements.b33 import B33

# sets
from sets.elementset import ElementSet

# materials (class-based composition)
from materials import (
    Material,
    ElasticityIsotropic,
    Density,
)

# beam profiles
from profiles.rect_profile import RectProfile

# sections
from sections.beam_section import BeamSection

# supports & loads
from supports.support import Support
from supports.support_collector import SupportCollector
from loads import CLoad, LoadCollector


def build_beam_model() -> Model:
    m = Model("Beam3Nodes")

    # --- NODES: 3 nodes along X ---
    n0 = m.add_node(Node(None, 0.0, 0.0, 0.0))
    n1 = m.add_node(Node(None, 1.0, 0.0, 0.0))
    n2 = m.add_node(Node(None, 2.0, 0.0, 0.0))

    # --- ELEMENTS: 2 beam elements (B33) ---
    e0 = m.add_element(B33(None, [n0, n1]))
    e1 = m.add_element(B33(None, [n1, n2]))

    # --- ELEMENT SET for section assignment ---
    beam_elset = m.add_elementset(ElementSet("BEAMS", [e0, e1]))

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
    m.add_section(BeamSection(profile=rect, material=steel, elset=beam_elset))

    # --- SUPPORTS ---
    # semantics: 0 = constrained, None = free (blank), 1 = prescribed "1"
    # left node fixed in all 6 DOFs
    sup_left = Support(n0, (0, 0, 0, 0, 0, 0))
    # middle node restrained only vertically (UY); others free
    sup_mid = Support(n1, (None, 0, None, None, None, None))

    sc = SupportCollector("BCS")
    sc.add(sup_left)
    sc.add(sup_mid)
    m.add_supportcollector(sc)

    # --- LOADS ---
    # concentrated load at the right node: downward 1 kN in Y
    lc = LoadCollector("LOADS", [CLoad(n2, (0.0, -1000.0, 0.0))])
    m.add_loadcollector(lc)

    return m


if __name__ == "__main__":
    model = build_beam_model()
    print(model.to_femaster())
