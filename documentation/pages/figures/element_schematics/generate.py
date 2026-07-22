from __future__ import annotations

"""
generate_element_schematics.py

Create clean schematic finite-element reference images for common FEM elements.
The script visualizes, in a documentation-friendly scientific style:

- node numbering in one large element view
- optional face numbering for solid elements
- integration point arrangement if applicable
- compact element properties

Supported element families:

- 3D solid elements: C3D4, C3D5, C3D6, C3D8, C3D10, C3D13, C3D15, C3D20, C3D20R
- shell elements: S3, S4, MITC4, MITC3FRT, MITC4FRT, S6, MITC6FRT, S8, MITC8, MITC8FRT, QSPT
- line elements: T3, B33

Dependencies:
    pip install numpy matplotlib pillow

Usage:
    python generate_element_schematics.py
    python generate_element_schematics.py --element C3D8
    python generate_element_schematics.py --element S4 --element B33
    python generate_element_schematics.py --all --out element_figures
    python generate_element_schematics.py --format svg
    python generate_element_schematics.py --show-faces
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Tuple
import argparse
import math
import textwrap

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from PIL import Image


Vec3 = Tuple[float, float, float]
Face = Tuple[int, ...]
Edge = Tuple[int, int]
ElementFamily = Literal["solid", "shell", "line"]


# -----------------------------------------------------------------------------
# Global visual style
# -----------------------------------------------------------------------------

rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "font.size": 11,
})

COLOR_EDGE = "0.08"
COLOR_EDGE_SECONDARY = "0.30"
COLOR_NODE = "0.02"
COLOR_MID_NODE = "0.14"
COLOR_FACE = (0.72, 0.72, 0.72, 0.13)
COLOR_FACE_EDGE = "0.55"
COLOR_IP = "0.05"
COLOR_LABEL_EDGE = "0.55"
COLOR_LABEL_FACE = "white"


@dataclass(frozen=True)
class ElementDefinition:
    name: str
    title: str
    family: ElementFamily
    nodes: Dict[int, Vec3]
    faces: Dict[int, Face]
    integration_domain: str
    integration_order: str
    shape_functions: str
    dof_per_node: int
    line_edges: Optional[Tuple[Edge, ...]] = None
    has_faces_in_solver: bool = True

    @property
    def node_count(self) -> int:
        return len(self.nodes)

    @property
    def face_count(self) -> int:
        if self.family != "solid" or not self.has_faces_in_solver:
            return 0
        return len(self.faces)

    @property
    def edge_count(self) -> int:
        if self.family == "line":
            return len(self.plot_edges)
        return len(edges_from_faces(self.faces.values()))

    @property
    def dof_count(self) -> int:
        return self.node_count * self.dof_per_node

    @property
    def plot_edges(self) -> List[Edge]:
        if self.line_edges is not None:
            return list(self.line_edges)
        if self.family in {"solid", "shell"}:
            return edges_from_faces(self.faces.values())
        return []


# -----------------------------------------------------------------------------
# Element definitions
# -----------------------------------------------------------------------------

ELEMENTS: Dict[str, ElementDefinition] = {
    "C3D4": ElementDefinition(
        name="C3D4",
        title="4-node tetrahedral solid element",
        family="solid",
        nodes={
            1: (1.0, 0.0, 0.0),
            2: (0.0, 1.0, 0.0),
            3: (0.0, 0.0, 0.0),
            4: (0.0, 0.0, 1.0),
        },
        faces={
            1: (1, 2, 3),
            2: (1, 4, 2),
            3: (2, 4, 3),
            4: (3, 4, 1),
        },
        integration_domain="DOMAIN_ISO_TET",
        integration_order="ORDER_LINEAR",
        shape_functions="linear tetrahedral Lagrange",
        dof_per_node=3,
    ),
    "C3D5": ElementDefinition(
        name="C3D5",
        title="5-node pyramid input element",
        family="solid",
        nodes={
            1: (1.0, 1.0, 0.0),
            2: (-1.0, 1.0, 0.0),
            3: (-1.0, -1.0, 0.0),
            4: (1.0, -1.0, 0.0),
            5: (0.0, 0.0, 1.15),
        },
        faces={
            1: (1, 2, 3, 4),
            2: (1, 5, 2),
            3: (2, 5, 3),
            4: (3, 5, 4),
            5: (4, 5, 1),
        },
        integration_domain="degenerated DOMAIN_ISO_HEX",
        integration_order="mapped C3D8 rule",
        shape_functions="pyramid input mapped to C3D8",
        dof_per_node=3,
        has_faces_in_solver=False,
    ),
    "C3D6": ElementDefinition(
        name="C3D6",
        title="6-node wedge solid element",
        family="solid",
        nodes={
            1: (1.0, 0.0, -1.0),
            2: (0.0, 1.0, -1.0),
            3: (0.0, 0.0, -1.0),
            4: (1.0, 0.0, 1.0),
            5: (0.0, 1.0, 1.0),
            6: (0.0, 0.0, 1.0),
        },
        faces={
            1: (1, 2, 3),
            2: (4, 5, 6),
            3: (1, 2, 5, 4),
            4: (2, 3, 6, 5),
            5: (3, 1, 4, 6),
        },
        integration_domain="DOMAIN_ISO_WEDGE",
        integration_order="ORDER_SUPER_LINEAR",
        shape_functions="linear wedge / prism",
        dof_per_node=3,
    ),
    "C3D8": ElementDefinition(
        name="C3D8",
        title="8-node hexahedral solid element",
        family="solid",
        nodes={
            1: (-1.0, -1.0, -1.0),
            2: (1.0, -1.0, -1.0),
            3: (1.0, 1.0, -1.0),
            4: (-1.0, 1.0, -1.0),
            5: (-1.0, -1.0, 1.0),
            6: (1.0, -1.0, 1.0),
            7: (1.0, 1.0, 1.0),
            8: (-1.0, 1.0, 1.0),
        },
        faces={
            1: (1, 2, 3, 4),
            2: (5, 6, 7, 8),
            3: (1, 2, 6, 5),
            4: (2, 3, 7, 6),
            5: (3, 4, 8, 7),
            6: (4, 1, 5, 8),
        },
        integration_domain="DOMAIN_ISO_HEX",
        integration_order="ORDER_QUADRATIC",
        shape_functions="trilinear Lagrange",
        dof_per_node=3,
    ),
    "C3D8R": ElementDefinition(
        name="C3D8R",
        title="8-node reduced-integration hexahedral solid element",
        family="solid",
        nodes={
            1: (-1.0, -1.0, -1.0),
            2: (1.0, -1.0, -1.0),
            3: (1.0, 1.0, -1.0),
            4: (-1.0, 1.0, -1.0),
            5: (-1.0, -1.0, 1.0),
            6: (1.0, -1.0, 1.0),
            7: (1.0, 1.0, 1.0),
            8: (-1.0, 1.0, 1.0),
        },
        faces={
            1: (1, 2, 3, 4),
            2: (5, 6, 7, 8),
            3: (1, 2, 6, 5),
            4: (2, 3, 7, 6),
            5: (3, 4, 8, 7),
            6: (4, 1, 5, 8),
        },
        integration_domain="DOMAIN_ISO_HEX",
        integration_order="ORDER_CONSTANT",
        shape_functions="trilinear Lagrange with hourglass stabilization",
        dof_per_node=3,
    ),
    "C3D10": ElementDefinition(
        name="C3D10",
        title="10-node tetrahedral solid element",
        family="solid",
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.0, 1.0, 0.0),
            4: (0.0, 0.0, 1.0),
            5: (0.5, 0.0, 0.0),
            6: (0.5, 0.5, 0.0),
            7: (0.0, 0.5, 0.0),
            8: (0.0, 0.0, 0.5),
            9: (0.5, 0.0, 0.5),
            10: (0.0, 0.5, 0.5),
        },
        faces={
            1: (1, 2, 3, 5, 6, 7),
            2: (1, 4, 2, 8, 9, 5),
            3: (2, 4, 3, 9, 10, 6),
            4: (3, 4, 1, 10, 8, 7),
        },
        integration_domain="DOMAIN_ISO_TET",
        integration_order="ORDER_QUADRATIC",
        shape_functions="quadratic tetrahedral Lagrange",
        dof_per_node=3,
    ),
    "C3D13": ElementDefinition(
        name="C3D13",
        title="13-node pyramid solid element",
        family="solid",
        nodes={
            1: (1.0, 1.0, 0.0),
            2: (-1.0, 1.0, 0.0),
            3: (-1.0, -1.0, 0.0),
            4: (1.0, -1.0, 0.0),
            5: (0.0, 0.0, 1.0),
            6: (0.0, 1.0, 0.0),
            7: (-1.0, 0.0, 0.0),
            8: (0.0, -1.0, 0.0),
            9: (1.0, 0.0, 0.0),
            10: (0.5, 0.5, 0.5),
            11: (-0.5, 0.5, 0.5),
            12: (-0.5, -0.5, 0.5),
            13: (0.5, -0.5, 0.5),
        },
        # These faces are used only for drawing edges. The supplied solver code does
        # not expose implemented surface definitions for C3D13.
        faces={
            1: (1, 2, 3, 4, 6, 7, 8, 9),
            2: (1, 5, 2, 10, 11, 6),
            3: (2, 5, 3, 11, 12, 7),
            4: (3, 5, 4, 12, 13, 8),
            5: (4, 5, 1, 13, 10, 9),
        },
        integration_domain="DOMAIN_ISO_PYRAMID",
        integration_order="ORDER_QUARTIC",
        shape_functions="quadratic pyramid",
        dof_per_node=3,
        has_faces_in_solver=False,
    ),
    "C3D15": ElementDefinition(
        name="C3D15",
        title="15-node wedge solid element",
        family="solid",
        nodes={
            1: (0.0, 0.0, -1.0),
            2: (1.0, 0.0, -1.0),
            3: (0.0, 1.0, -1.0),
            4: (0.0, 0.0, 1.0),
            5: (1.0, 0.0, 1.0),
            6: (0.0, 1.0, 1.0),
            7: (0.5, 0.0, -1.0),
            8: (0.5, 0.5, -1.0),
            9: (0.0, 0.5, -1.0),
            10: (0.5, 0.0, 1.0),
            11: (0.5, 0.5, 1.0),
            12: (0.0, 0.5, 1.0),
            13: (0.0, 0.0, 0.0),
            14: (1.0, 0.0, 0.0),
            15: (0.0, 1.0, 0.0),
        },
        faces={
            1: (1, 2, 3, 7, 8, 9),
            2: (4, 5, 6, 10, 11, 12),
            3: (1, 2, 5, 4, 7, 14, 10, 13),
            4: (2, 3, 6, 5, 8, 15, 11, 14),
            5: (3, 1, 4, 6, 9, 13, 12, 15),
        },
        integration_domain="DOMAIN_ISO_WEDGE",
        integration_order="ORDER_SUPER_QUADRATIC",
        shape_functions="quadratic wedge / prism",
        dof_per_node=3,
    ),
    "C3D20": ElementDefinition(
        name="C3D20",
        title="20-node hexahedral solid element",
        family="solid",
        nodes={
            1: (-1.0, -1.0, -1.0),
            2: (1.0, -1.0, -1.0),
            3: (1.0, 1.0, -1.0),
            4: (-1.0, 1.0, -1.0),
            5: (-1.0, -1.0, 1.0),
            6: (1.0, -1.0, 1.0),
            7: (1.0, 1.0, 1.0),
            8: (-1.0, 1.0, 1.0),
            9: (0.0, -1.0, -1.0),
            10: (1.0, 0.0, -1.0),
            11: (0.0, 1.0, -1.0),
            12: (-1.0, 0.0, -1.0),
            13: (0.0, -1.0, 1.0),
            14: (1.0, 0.0, 1.0),
            15: (0.0, 1.0, 1.0),
            16: (-1.0, 0.0, 1.0),
            17: (-1.0, -1.0, 0.0),
            18: (1.0, -1.0, 0.0),
            19: (1.0, 1.0, 0.0),
            20: (-1.0, 1.0, 0.0),
        },
        faces={
            1: (1, 2, 3, 4, 9, 10, 11, 12),
            2: (5, 6, 7, 8, 13, 14, 15, 16),
            3: (1, 2, 6, 5, 9, 18, 13, 17),
            4: (2, 3, 7, 6, 10, 19, 14, 18),
            5: (3, 4, 8, 7, 11, 20, 15, 19),
            6: (4, 1, 5, 8, 12, 17, 16, 20),
        },
        integration_domain="DOMAIN_ISO_HEX",
        integration_order="ORDER_QUARTIC",
        shape_functions="quadratic serendipity hexahedral",
        dof_per_node=3,
    ),
    "C3D20R": ElementDefinition(
        name="C3D20R",
        title="20-node reduced-integration hexahedral solid element",
        family="solid",
        nodes={
            1: (-1.0, -1.0, -1.0),
            2: (1.0, -1.0, -1.0),
            3: (1.0, 1.0, -1.0),
            4: (-1.0, 1.0, -1.0),
            5: (-1.0, -1.0, 1.0),
            6: (1.0, -1.0, 1.0),
            7: (1.0, 1.0, 1.0),
            8: (-1.0, 1.0, 1.0),
            9: (0.0, -1.0, -1.0),
            10: (1.0, 0.0, -1.0),
            11: (0.0, 1.0, -1.0),
            12: (-1.0, 0.0, -1.0),
            13: (0.0, -1.0, 1.0),
            14: (1.0, 0.0, 1.0),
            15: (0.0, 1.0, 1.0),
            16: (-1.0, 0.0, 1.0),
            17: (-1.0, -1.0, 0.0),
            18: (1.0, -1.0, 0.0),
            19: (1.0, 1.0, 0.0),
            20: (-1.0, 1.0, 0.0),
        },
        faces={
            1: (1, 2, 3, 4, 9, 10, 11, 12),
            2: (5, 6, 7, 8, 13, 14, 15, 16),
            3: (1, 2, 6, 5, 9, 18, 13, 17),
            4: (2, 3, 7, 6, 10, 19, 14, 18),
            5: (3, 4, 8, 7, 11, 20, 15, 19),
            6: (4, 1, 5, 8, 12, 17, 16, 20),
        },
        integration_domain="DOMAIN_ISO_HEX",
        integration_order="ORDER_QUADRATIC",
        shape_functions="quadratic serendipity hexahedral",
        dof_per_node=3,
    ),
    "S3": ElementDefinition(
        name="S3",
        title="3-node triangular shell element",
        family="shell",
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.20, 0.88, 0.0),
        },
        faces={1: (1, 2, 3)},
        integration_domain="DOMAIN_ISO_TRI",
        integration_order="ORDER_CUBIC",
        shape_functions="linear triangular shell",
        dof_per_node=6,
    ),
    "S4": ElementDefinition(
        name="S4",
        title="4-node quadrilateral shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
        },
        faces={1: (1, 2, 3, 4)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_CUBIC",
        shape_functions="bilinear quadrilateral shell",
        dof_per_node=6,
    ),
    "MITC4": ElementDefinition(
        name="MITC4",
        title="4-node MITC quadrilateral shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
        },
        faces={1: (1, 2, 3, 4)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_CUBIC",
        shape_functions="bilinear MITC quadrilateral shell",
        dof_per_node=6,
    ),
    "MITC3FRT": ElementDefinition(
        name="MITC3FRT",
        title="3-node finite-rotation MITC shell element",
        family="shell",
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.20, 0.88, 0.0),
        },
        faces={1: (1, 2, 3)},
        integration_domain="DOMAIN_ISO_TRI",
        integration_order="ORDER_CUBIC",
        shape_functions="finite-rotation MITC3 shell",
        dof_per_node=6,
    ),
    "MITC4FRT": ElementDefinition(
        name="MITC4FRT",
        title="4-node finite-rotation MITC shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
        },
        faces={1: (1, 2, 3, 4)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_CUBIC",
        shape_functions="finite-rotation MITC4 shell",
        dof_per_node=6,
    ),
    "S6": ElementDefinition(
        name="S6",
        title="6-node triangular shell element",
        family="shell",
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.20, 0.88, 0.0),
            4: (0.5, 0.0, 0.0),
            5: (0.60, 0.44, 0.0),
            6: (0.10, 0.44, 0.0),
        },
        faces={1: (1, 2, 3, 4, 5, 6)},
        integration_domain="DOMAIN_ISO_TRI",
        integration_order="ORDER_CUBIC",
        shape_functions="quadratic triangular shell",
        dof_per_node=6,
    ),
    "MITC6FRT": ElementDefinition(
        name="MITC6FRT",
        title="6-node finite-rotation MITC6-b shell element",
        family="shell",
        nodes={
            1: (0.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
            3: (0.20, 0.88, 0.0),
            4: (0.5, 0.0, 0.0),
            5: (0.60, 0.44, 0.0),
            6: (0.10, 0.44, 0.0),
        },
        faces={1: (1, 2, 3, 4, 5, 6)},
        integration_domain="DOMAIN_ISO_TRI",
        integration_order="ORDER_CUBIC",
        shape_functions="finite-rotation MITC6-b shell",
        dof_per_node=6,
    ),
    "S8": ElementDefinition(
        name="S8",
        title="8-node quadrilateral shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
            5: (0.0, -1.0, 0.0),
            6: (1.0, 0.0, 0.0),
            7: (0.0, 1.0, 0.0),
            8: (-1.0, 0.0, 0.0),
        },
        faces={1: (1, 2, 3, 4, 5, 6, 7, 8)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_QUINTIC",
        shape_functions="quadratic serendipity shell",
        dof_per_node=6,
    ),
    "MITC8": ElementDefinition(
        name="MITC8",
        title="8-node MITC quadrilateral shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
            5: (0.0, -1.0, 0.0),
            6: (1.0, 0.0, 0.0),
            7: (0.0, 1.0, 0.0),
            8: (-1.0, 0.0, 0.0),
        },
        faces={1: (1, 2, 3, 4, 5, 6, 7, 8)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_QUINTIC",
        shape_functions="quadratic MITC8 shell",
        dof_per_node=6,
    ),
    "MITC8FRT": ElementDefinition(
        name="MITC8FRT",
        title="8-node finite-rotation MITC shell element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
            5: (0.0, -1.0, 0.0),
            6: (1.0, 0.0, 0.0),
            7: (0.0, 1.0, 0.0),
            8: (-1.0, 0.0, 0.0),
        },
        faces={1: (1, 2, 3, 4, 5, 6, 7, 8)},
        integration_domain="DOMAIN_ISO_QUAD",
        integration_order="ORDER_QUINTIC",
        shape_functions="finite-rotation MITC8 shell",
        dof_per_node=6,
    ),
    "QSPT": ElementDefinition(
        name="QSPT",
        title="4-node quadrilateral shear-panel element",
        family="shell",
        nodes={
            1: (-1.0, -1.0, 0.0),
            2: (1.0, -1.0, 0.0),
            3: (1.0, 1.0, 0.0),
            4: (-1.0, 1.0, 0.0),
        },
        faces={1: (1, 2, 3, 4)},
        integration_domain="DOMAIN_NONE",
        integration_order="ORDER_NONE",
        shape_functions="quadrilateral shear-panel formulation",
        dof_per_node=3,
    ),
    "T3": ElementDefinition(
        name="T3",
        title="2-node truss element",
        family="line",
        nodes={
            1: (-1.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
        },
        faces={},
        line_edges=((1, 2),),
        integration_domain="DOMAIN_ISO_LINE",
        integration_order="ORDER_LINEAR",
        shape_functions="linear 2-node truss",
        dof_per_node=3,
    ),
    "B33": ElementDefinition(
        name="B33",
        title="2-node beam element",
        family="line",
        nodes={
            1: (-1.0, 0.0, 0.0),
            2: (1.0, 0.0, 0.0),
        },
        faces={},
        line_edges=((1, 2),),
        integration_domain="DOMAIN_ISO_LINE",
        integration_order="ORDER_QUADRATIC",
        shape_functions="linear 2-node beam",
        dof_per_node=6,
    ),
}


# -----------------------------------------------------------------------------
# Connectivity helpers
# -----------------------------------------------------------------------------

def face_corner_nodes(face: Face) -> Face:
    """Return geometric corner nodes from a face definition."""
    if len(face) in (2, 3, 4):
        return face
    if len(face) == 6:
        return face[:3]
    if len(face) == 8:
        return face[:4]
    raise ValueError(f"Unsupported face length {len(face)} for face {face}")


def face_midside_edges(face: Face) -> List[Tuple[int, int, int]]:
    """Return corner-midside-corner triples for quadratic faces."""
    if len(face) == 6:
        c1, c2, c3, m12, m23, m31 = face
        return [(c1, m12, c2), (c2, m23, c3), (c3, m31, c1)]
    if len(face) == 8:
        c1, c2, c3, c4, m12, m23, m34, m41 = face
        return [(c1, m12, c2), (c2, m23, c3), (c3, m34, c4), (c4, m41, c1)]
    return []


def edges_from_faces(faces: Iterable[Face]) -> List[Edge]:
    edges: set[Edge] = set()
    for face in faces:
        corners = face_corner_nodes(face)
        for a, b in zip(corners, corners[1:] + corners[:1]):
            edges.add(tuple(sorted((a, b))))
    return sorted(edges)


# -----------------------------------------------------------------------------
# Quadrature points for visual documentation
# -----------------------------------------------------------------------------

def gauss_legendre_1d(order: str) -> List[float]:
    if order in {"ORDER_CONSTANT", "ORDER_LINEAR"}:
        return [0.0]
    if order in {"ORDER_QUADRATIC", "ORDER_CUBIC", "ORDER_SUPER_LINEAR"}:
        return [-1.0 / math.sqrt(3.0), 1.0 / math.sqrt(3.0)]
    if order in {"ORDER_QUARTIC", "ORDER_QUINTIC", "ORDER_SUPER_QUADRATIC"}:
        return [-math.sqrt(3.0 / 5.0), 0.0, math.sqrt(3.0 / 5.0)]
    return []


def integration_points(defn: ElementDefinition) -> np.ndarray:
    """Return integration points in reference coordinates.

    The point order follows the FEMaster quadrature registrations:

    - DOMAIN_ISO_QUAD:
      ORDER_CONSTANT/LINEAR  -> 1 point
      ORDER_QUADRATIC/CUBIC  -> 4 points
      ORDER_QUARTIC/QUINTIC  -> 9 points
    - DOMAIN_ISO_TRI:
      ORDER_CONSTANT/LINEAR  -> 1 point
      ORDER_QUADRATIC        -> 3 points
      ORDER_CUBIC            -> 4 points

    The third coordinate is used as z/t for 3D domains and remains zero for
    shell/line reference domains.
    """
    domain = defn.integration_domain
    order = defn.integration_order

    if domain == "DOMAIN_NONE" or order == "ORDER_NONE":
        return np.zeros((0, 3), dtype=float)

    if domain == "DOMAIN_ISO_LINE":
        pts = gauss_legendre_1d(order)
        return np.array([(r, 0.0, 0.0) for r in pts], dtype=float)

    if domain == "DOMAIN_ISO_QUAD":
        if order in {"ORDER_CONSTANT", "ORDER_LINEAR"}:
            return np.array([[0.0, 0.0, 0.0]], dtype=float)
        if order in {"ORDER_QUADRATIC", "ORDER_CUBIC"}:
            a = math.sqrt(3.0) / 3.0
            return np.array([
                [ a,  a, 0.0],
                [-a,  a, 0.0],
                [-a, -a, 0.0],
                [ a, -a, 0.0],
            ], dtype=float)
        if order in {"ORDER_QUARTIC", "ORDER_QUINTIC"}:
            a = math.sqrt(0.6)
            return np.array([
                [-a, -a, 0.0], [0.0, -a, 0.0], [a, -a, 0.0],
                [-a, 0.0, 0.0], [0.0, 0.0, 0.0], [a, 0.0, 0.0],
                [-a, a, 0.0], [0.0, a, 0.0], [a, a, 0.0],
            ], dtype=float)
        return np.zeros((0, 3), dtype=float)

    if domain == "DOMAIN_ISO_TRI":
        if order in {"ORDER_CONSTANT", "ORDER_LINEAR"}:
            return np.array([[1.0 / 3.0, 1.0 / 3.0, 0.0]], dtype=float)
        if order == "ORDER_QUADRATIC":
            return np.array([
                [1.0 / 6.0, 1.0 / 6.0, 0.0],
                [2.0 / 3.0, 1.0 / 6.0, 0.0],
                [1.0 / 6.0, 2.0 / 3.0, 0.0],
            ], dtype=float)
        if order == "ORDER_CUBIC":
            return np.array([
                [1.0 / 3.0, 1.0 / 3.0, 0.0],
                [1.0 / 5.0, 3.0 / 5.0, 0.0],
                [1.0 / 5.0, 1.0 / 5.0, 0.0],
                [3.0 / 5.0, 1.0 / 5.0, 0.0],
            ], dtype=float)
        return np.zeros((0, 3), dtype=float)

    if domain == "DOMAIN_ISO_HEX":
        pts = gauss_legendre_1d(order)
        return np.array([(r, s, t) for r in pts for s in pts for t in pts], dtype=float)

    if domain == "DOMAIN_ISO_TET":
        if order in {"ORDER_CONSTANT", "ORDER_LINEAR"}:
            return np.array([[0.25, 0.25, 0.25]], dtype=float)
        a = 0.5854101966249685
        b = 0.1381966011250105
        return np.array([[a, b, b], [b, a, b], [b, b, a], [b, b, b]], dtype=float)

    if domain == "DOMAIN_ISO_WEDGE":
        if order == "ORDER_SUPER_LINEAR":
            tri = np.array([[1.0 / 3.0, 1.0 / 3.0]], dtype=float)
        else:
            tri = np.array([
                [1.0 / 6.0, 1.0 / 6.0],
                [2.0 / 3.0, 1.0 / 6.0],
                [1.0 / 6.0, 2.0 / 3.0],
            ], dtype=float)
        z = gauss_legendre_1d(order)
        return np.array([(r, s, t) for r, s in tri for t in z], dtype=float)

    if domain == "DOMAIN_ISO_PYRAMID":
        a = 0.45
        return np.array([
            [-a, -a, 0.25], [a, -a, 0.25], [a, a, 0.25], [-a, a, 0.25],
            [-a / 2.0, -a / 2.0, 0.65], [a / 2.0, -a / 2.0, 0.65],
            [a / 2.0, a / 2.0, 0.65], [-a / 2.0, a / 2.0, 0.65],
        ], dtype=float)

    return np.zeros((0, 3), dtype=float)


# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------

def set_equal_3d(ax, pts: np.ndarray, pad: float = 0.18, zoom: float = 1.0) -> None:
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = 0.5 * (mins + maxs)
    radius = 0.5 * float(np.max(maxs - mins)) + pad
    radius /= zoom

    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))
    ax.margins(0.0)


def element_zoom(defn: ElementDefinition) -> float:
    if defn.family == "line":
        return 1.20
    if defn.family == "shell":
        return 1.10
    if defn.name in {"C3D4", "C3D10"}:
        return 1.08
    if defn.name in {"C3D6", "C3D15", "C3D13"}:
        return 1.06
    return 1.04


def style_3d_axis(ax, defn: ElementDefinition) -> None:
    ax.set_axis_off()
    if defn.family == "line":
        ax.view_init(elev=16, azim=-62)
    elif defn.family == "shell":
        ax.view_init(elev=28, azim=-56)
    else:
        ax.view_init(elev=22, azim=-53)


def label_offset(p: np.ndarray, all_points: np.ndarray, scale: float = 0.11) -> np.ndarray:
    center = all_points.mean(axis=0)
    v = p - center
    n = np.linalg.norm(v)
    if n < 1e-12:
        return np.array([scale, scale, scale])
    return scale * v / n


def transformed_nodes(defn: ElementDefinition) -> Dict[int, Vec3]:
    """Return coordinates used only for plotting."""
    z_scale = 1.45 if defn.name == "C3D13" else 1.0
    return {k: (x, y, z * z_scale) for k, (x, y, z) in defn.nodes.items()}


def shell_shape_position(defn: ElementDefinition, rst: np.ndarray) -> np.ndarray:
    """Map shell reference coordinates to the plotted element geometry."""
    r = float(rst[0])
    s = float(rst[1])
    nodes = transformed_nodes(defn)

    if defn.integration_domain == "DOMAIN_ISO_TRI":
        p1 = np.array(nodes[1], dtype=float)
        p2 = np.array(nodes[2], dtype=float)
        p3 = np.array(nodes[3], dtype=float)
        if defn.node_count == 3:
            return (1.0 - r - s) * p1 + r * p2 + s * p3

        # Quadratic triangular Lagrange interpolation. The local convention is
        # N1=(1-r-s)(1-2r-2s), N2=r(2r-1), N3=s(2s-1),
        # N4=4r(1-r-s), N5=4rs, N6=4s(1-r-s).
        l1 = 1.0 - r - s
        l2 = r
        l3 = s
        shape = np.array([
            l1 * (2.0 * l1 - 1.0),
            l2 * (2.0 * l2 - 1.0),
            l3 * (2.0 * l3 - 1.0),
            4.0 * l1 * l2,
            4.0 * l2 * l3,
            4.0 * l3 * l1,
            ])
        coords = np.array([nodes[i] for i in range(1, 7)], dtype=float)
        return shape @ coords

    if defn.integration_domain == "DOMAIN_ISO_QUAD":
        p = {i: np.array(nodes[i], dtype=float) for i in nodes}
        if defn.node_count == 4:
            shape = {
                1: 0.25 * (1.0 - r) * (1.0 - s),
                2: 0.25 * (1.0 + r) * (1.0 - s),
                3: 0.25 * (1.0 + r) * (1.0 + s),
                4: 0.25 * (1.0 - r) * (1.0 + s),
            }
            return sum(shape[i] * p[i] for i in range(1, 5))

        # 8-node serendipity interpolation.
        shape = {
            1: -0.25 * (1.0 - r) * (1.0 - s) * (1.0 + r + s),
            2: -0.25 * (1.0 + r) * (1.0 - s) * (1.0 - r + s),
            3: -0.25 * (1.0 + r) * (1.0 + s) * (1.0 - r - s),
            4: -0.25 * (1.0 - r) * (1.0 + s) * (1.0 + r - s),
            5: 0.5 * (1.0 - r * r) * (1.0 - s),
            6: 0.5 * (1.0 + r) * (1.0 - s * s),
            7: 0.5 * (1.0 - r * r) * (1.0 + s),
            8: 0.5 * (1.0 - r) * (1.0 - s * s),
        }
        return sum(shape[i] * p[i] for i in range(1, min(defn.node_count, 8) + 1))

    return rst.copy()


def transformed_integration_points(defn: ElementDefinition) -> np.ndarray:
    ips = integration_points(defn).copy()
    if len(ips) == 0:
        return ips

    if defn.family == "shell":
        return np.array([shell_shape_position(defn, ip) for ip in ips], dtype=float)

    z_scale = 1.45 if defn.name == "C3D13" else 1.0
    ips[:, 2] *= z_scale
    return ips


def node_marker_size(defn: ElementDefinition, node_id: int) -> float:
    if defn.family == "line":
        return 58.0
    if defn.family == "shell" and node_id > len(face_corner_nodes(next(iter(defn.faces.values())))):
        return 38.0
    if defn.family == "solid" and node_id > 8:
        return 38.0
    return 48.0


def draw_curved_quadratic_edges(ax, defn: ElementDefinition, nodes: Dict[int, Vec3]) -> None:
    """Draw quadratic element midside connectivity with subtle segmented edges."""
    segments = []
    seen = set()
    for face in defn.faces.values():
        for a, m, b in face_midside_edges(face):
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            pa = np.array(nodes[a], dtype=float)
            pm = np.array(nodes[m], dtype=float)
            pb = np.array(nodes[b], dtype=float)
            segments.append([pa, pm])
            segments.append([pm, pb])

    if not segments:
        return

    collection = Line3DCollection(segments, colors=COLOR_EDGE, linewidths=1.45)
    ax.add_collection3d(collection)


def draw_element(
        ax,
        defn: ElementDefinition,
        *,
        show_nodes: bool = True,
        show_node_labels: bool = True,
        show_faces: bool = True,
        show_face_labels: bool = False,
        show_integration_points: bool = True,
        show_face_leaders: bool = True,
) -> None:
    nodes = transformed_nodes(defn)
    coords = np.array([nodes[i] for i in sorted(nodes)], dtype=float)

    if defn.family in {"solid", "shell"} and show_faces:
        polys = []
        for face in defn.faces.values():
            corners = face_corner_nodes(face)
            polys.append([nodes[i] for i in corners])
        collection = Poly3DCollection(polys, alpha=0.13, linewidths=0.65, edgecolors=COLOR_FACE_EDGE)
        collection.set_facecolor(COLOR_FACE)
        ax.add_collection3d(collection)

    # Draw primary/corner edges.
    for a, b in defn.plot_edges:
        pa = np.array(nodes[a])
        pb = np.array(nodes[b])
        ax.plot(
            [pa[0], pb[0]],
            [pa[1], pb[1]],
            [pa[2], pb[2]],
            color=COLOR_EDGE,
            lw=1.55 if defn.family != "line" else 2.05,
            solid_capstyle="round",
        )

    # Draw midside connectivity for quadratic shells/solids.
    if defn.family in {"solid", "shell"}:
        draw_curved_quadratic_edges(ax, defn, nodes)

    # Beam/truss centerline direction marker.
    if defn.family == "line":
        pa = np.array(nodes[1], dtype=float)
        pb = np.array(nodes[2], dtype=float)
        mid = 0.5 * (pa + pb)
        direction = pb - pa
        direction /= np.linalg.norm(direction)
        ax.quiver(
            mid[0], mid[1], mid[2],
            0.34 * direction[0], 0.34 * direction[1], 0.34 * direction[2],
            arrow_length_ratio=0.24,
            color=COLOR_EDGE_SECONDARY,
            linewidth=1.1,
            )

    if show_integration_points:
        ips = transformed_integration_points(defn)
        if len(ips) > 0:
            ax.scatter(
                ips[:, 0], ips[:, 1], ips[:, 2],
                s=28 if defn.family != "line" else 36,
                color=COLOR_IP,
                marker="x",
                linewidths=1.2,
                depthshade=False,
                zorder=8,
            )

    if show_nodes:
        for node_id, p in nodes.items():
            p_arr = np.array(p, dtype=float)
            ax.scatter(
                [p_arr[0]], [p_arr[1]], [p_arr[2]],
                s=node_marker_size(defn, node_id),
                color=COLOR_NODE if node_id <= 8 else COLOR_MID_NODE,
                depthshade=False,
                zorder=9,
            )

    if show_node_labels:
        for i, p in nodes.items():
            p_arr = np.array(p, dtype=float)
            label_pos = p_arr + label_offset(p_arr, coords, 0.16 if defn.family != "line" else 0.13)
            ax.text(
                label_pos[0], label_pos[1], label_pos[2],
                str(i),
                fontsize=10.5,
                color="black",
                ha="center",
                va="center",
                bbox=dict(
                    boxstyle="square,pad=0.13",
                    fc=COLOR_LABEL_FACE,
                    ec=COLOR_LABEL_EDGE,
                    lw=0.45,
                    alpha=0.96,
                ),
            )

    if show_face_labels and defn.family == "solid" and defn.has_faces_in_solver:
        element_center = coords.mean(axis=0)
        for face_id, face in defn.faces.items():
            corners = np.array([nodes[i] for i in face_corner_nodes(face)], dtype=float)
            face_center = corners.mean(axis=0)
            outward = face_center - element_center
            nrm = np.linalg.norm(outward)
            if nrm < 1e-12:
                outward = np.array([0.0, 0.0, 1.0])
                nrm = 1.0
            outward = outward / nrm
            label_pos = face_center + 0.27 * outward

            if show_face_leaders:
                ax.plot(
                    [face_center[0], label_pos[0]],
                    [face_center[1], label_pos[1]],
                    [face_center[2], label_pos[2]],
                    color="0.42",
                    lw=0.75,
                    ls="-",
                    zorder=7,
                )

            ax.text(
                label_pos[0], label_pos[1], label_pos[2],
                f"F{face_id}",
                fontsize=10.0,
                ha="center",
                va="center",
                color="0.10",
                bbox=dict(
                    boxstyle="square,pad=0.14",
                    fc="white",
                    ec="0.45",
                    lw=0.45,
                    alpha=0.95,
                ),
            )

    set_equal_3d(ax, coords, pad=0.22, zoom=element_zoom(defn))
    style_3d_axis(ax, defn)


def crop_rgba_to_content(image: np.ndarray, alpha_threshold: int = 8, padding: int = 12) -> np.ndarray:
    """Crop an RGBA image to its visible, non-transparent content."""
    if image.ndim != 3 or image.shape[2] != 4:
        raise ValueError("Expected an RGBA image with shape (H, W, 4).")

    alpha = image[:, :, 3]
    ys, xs = np.where(alpha > alpha_threshold)
    if len(xs) == 0 or len(ys) == 0:
        return image

    h, w = alpha.shape
    x0 = max(int(xs.min()) - padding, 0)
    x1 = min(int(xs.max()) + padding + 1, w)
    y0 = max(int(ys.min()) - padding, 0)
    y1 = min(int(ys.max()) + padding + 1, h)
    return image[y0:y1, x0:x1, :]


def render_element_rgba(
        defn: ElementDefinition,
        *,
        dpi: int = 300,
        show_faces: bool = True,
        show_face_labels: bool = False,
        show_integration_points: bool = True,
) -> np.ndarray:
    """Render only the element to a transparent RGBA image and crop it."""
    fig = plt.figure(figsize=(7.2, 7.2), dpi=dpi)
    fig.patch.set_alpha(0.0)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection="3d")
    ax.patch.set_alpha(0.0)

    draw_element(
        ax,
        defn,
        show_nodes=True,
        show_node_labels=True,
        show_faces=show_faces,
        show_face_labels=show_face_labels,
        show_face_leaders=show_face_labels,
        show_integration_points=show_integration_points,
    )

    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    buffer = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    image = buffer.reshape((height, width, 4)).copy()
    plt.close(fig)

    return crop_rgba_to_content(image, padding=20)


def write_properties(ax, defn: ElementDefinition) -> None:
    ax.axis("off")

    family_name = {
        "solid": "solid",
        "shell": "shell",
        "line": "line",
    }[defn.family]

    rows = [
        ("Element", defn.name),
        ("Type", family_name),
        ("Nodes", str(defn.node_count)),
        ("Faces", "-" if defn.face_count == 0 else str(defn.face_count)),
        ("Edges", str(defn.edge_count)),
        ("DOF / node", str(defn.dof_per_node)),
        ("Total DOF", str(defn.dof_count)),
        ("Shape functions", "\n".join(textwrap.wrap(defn.shape_functions, width=28))),
        ("Integration points", str(len(integration_points(defn))) if len(integration_points(defn)) > 0 else "-"),
    ]

    table_height = 0.52 if defn.family != "line" else 0.48
    table_y0 = 0.5 - table_height / 2.0

    table = ax.table(
        cellText=[[k, v] for k, v in rows],
        cellLoc="left",
        colLoc="left",
        loc="center",
        bbox=[0.0, table_y0, 1.00, table_height],
        colWidths=[0.36, 0.64],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.8)
    table.scale(1.0, 1.18)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("0.70")
        cell.set_linewidth(0.45)
        cell.set_facecolor("white")
        if col == 0:
            cell.set_text_props(weight="normal", color="0.18")
        else:
            cell.set_text_props(color="0.04")
            cell.get_text().set_wrap(True)


def create_element_figure(
        defn: ElementDefinition,
        out_path: Path,
        *,
        dpi: int = 300,
        show_faces: bool = True,
        show_face_labels: bool = False,
        show_integration_points: bool = True,
) -> None:
    element_image = render_element_rgba(
        defn,
        dpi=dpi,
        show_faces=show_faces,
        show_face_labels=show_face_labels,
        show_integration_points=show_integration_points,
    )

    fig = plt.figure(figsize=(12.0, 6.2), constrained_layout=False)
    gs = fig.add_gridspec(
        1, 2,
        width_ratios=[1.55, 0.95],
        left=0.030,
        right=0.975,
        top=0.965,
        bottom=0.055,
        wspace=0.030,
    )

    # ------------------------------------------------------------------
    # Left: fixed-size element image area
    # ------------------------------------------------------------------

    ax_img = fig.add_subplot(gs[0, 0])
    ax_img.axis("off")

    canvas_h = 1500
    canvas_w = 1500
    canvas = np.zeros((canvas_h, canvas_w, 4), dtype=np.uint8)

    img_h, img_w = element_image.shape[:2]

    pad_top_bottom = 115
    pad_left_right = 95
    max_h = canvas_h - 2 * pad_top_bottom
    max_w = canvas_w - 2 * pad_left_right

    scale = min(max_w / img_w, max_h / img_h)
    new_w = int(img_w * scale)
    new_h = int(img_h * scale)

    resized = np.array(
        Image.fromarray(element_image).resize((new_w, new_h), Image.Resampling.LANCZOS)
    )

    y0 = (canvas_h - new_h) // 2
    x0 = (canvas_w - new_w) // 2
    canvas[y0:y0 + new_h, x0:x0 + new_w] = resized

    ax_img.imshow(canvas)
    ax_img.set_anchor("C")

    # ------------------------------------------------------------------
    # Right: vertically centered property table
    # ------------------------------------------------------------------

    ax_props = fig.add_subplot(gs[0, 1])
    ax_props.set_anchor("C")
    write_properties(ax_props, defn)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)


def generate(
        elements: Sequence[str],
        out_dir: Path,
        file_format: str,
        dpi: int,
        *,
        show_faces: bool,
        show_face_labels: bool,
        show_integration_points: bool,
) -> None:
    for name in elements:
        key = name.upper()
        if key not in ELEMENTS:
            valid = ", ".join(sorted(ELEMENTS))
            raise KeyError(f"Unknown element '{name}'. Valid elements: {valid}")
        defn = ELEMENTS[key]
        out = out_dir / f"{defn.name.lower()}_schematic.{file_format}"
        create_element_figure(
            defn,
            out,
            dpi=dpi,
            show_faces=show_faces,
            show_face_labels=show_face_labels,
            show_integration_points=show_integration_points,
        )
        print(f"wrote {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate schematic FEM element figures.")
    parser.add_argument(
        "--element",
        action="append",
        default=[],
        help="Element name to generate, e.g. C3D8, S4, B33. Can be passed multiple times.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate all available elements.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("element_schematics"),
        help="Output directory.",
    )
    parser.add_argument(
        "--format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Output image format.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for raster output.",
    )
    parser.add_argument(
        "--show-faces",
        action="store_true",
        help="Show translucent faces for shell/solid elements. Enabled by default for geometry clarity.",
    )
    parser.add_argument(
        "--hide-faces",
        action="store_true",
        help="Hide translucent element faces and draw only edges/nodes.",
    )
    parser.add_argument(
        "--show-face-labels",
        action="store_true",
        help="Show solid face labels F1, F2, ... . Disabled by default for a cleaner documentation style.",
    )
    parser.add_argument(
        "--hide-integration-points",
        action="store_true",
        help="Hide integration points.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = sorted(ELEMENTS) if args.all or not args.element else args.element

    show_faces = not args.hide_faces
    if args.show_faces:
        show_faces = True

    generate(
        selected,
        args.out,
        args.format,
        args.dpi,
        show_faces=show_faces,
        show_face_labels=args.show_face_labels,
        show_integration_points=not args.hide_integration_points,
    )


if __name__ == "__main__":
    main()
