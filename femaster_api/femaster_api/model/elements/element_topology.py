"""Element topology enum and FEMaster type mapping."""

from __future__ import annotations

from enum import Enum


class ElementTopology(Enum):
    TET4 = "TET4"
    PYRAMID5 = "PYRAMID5"
    WEDGE6 = "WEDGE6"
    HEX8 = "HEX8"
    TET10 = "TET10"
    WEDGE15 = "WEDGE15"
    HEX20 = "HEX20"
    HEX20R = "HEX20R"
    TRUSS2 = "TRUSS2"
    BEAM2 = "BEAM2"
    TRI3 = "TRI3"
    QUAD4 = "QUAD4"
    MITC4 = "MITC4"
    TRI6 = "TRI6"
    QUAD8 = "QUAD8"
    QSPT = "QSPT"


FEMASTER_ELEMENT_TYPES: dict[ElementTopology, str] = {
    ElementTopology.TET4: "C3D4",
    ElementTopology.PYRAMID5: "C3D5",
    ElementTopology.WEDGE6: "C3D6",
    ElementTopology.HEX8: "C3D8",
    ElementTopology.TET10: "C3D10",
    ElementTopology.WEDGE15: "C3D15",
    ElementTopology.HEX20: "C3D20",
    ElementTopology.HEX20R: "C3D20R",
    ElementTopology.TRUSS2: "T3",
    ElementTopology.BEAM2: "B33",
    ElementTopology.TRI3: "S3",
    ElementTopology.QUAD4: "S4",
    ElementTopology.MITC4: "MITC4",
    ElementTopology.TRI6: "S6",
    ElementTopology.QUAD8: "S8",
    ElementTopology.QSPT: "QSPT",
}

TOPOLOGY_NODE_COUNTS: dict[ElementTopology, tuple[int, ...]] = {
    ElementTopology.TET4: (4,),
    ElementTopology.PYRAMID5: (5,),
    ElementTopology.WEDGE6: (6,),
    ElementTopology.HEX8: (8,),
    ElementTopology.TET10: (10,),
    ElementTopology.WEDGE15: (15,),
    ElementTopology.HEX20: (20,),
    ElementTopology.HEX20R: (20,),
    ElementTopology.TRUSS2: (2,),
    ElementTopology.BEAM2: (2, 3),
    ElementTopology.TRI3: (3,),
    ElementTopology.QUAD4: (4,),
    ElementTopology.MITC4: (4,),
    ElementTopology.TRI6: (6,),
    ElementTopology.QUAD8: (8,),
    ElementTopology.QSPT: (4,),
}

C3D4 = ElementTopology.TET4
C3D5 = ElementTopology.PYRAMID5
C3D6 = ElementTopology.WEDGE6
C3D8 = ElementTopology.HEX8
C3D10 = ElementTopology.TET10
C3D15 = ElementTopology.WEDGE15
C3D20 = ElementTopology.HEX20
C3D20R = ElementTopology.HEX20R
T3 = ElementTopology.TRUSS2
B33 = ElementTopology.BEAM2
S3 = ElementTopology.TRI3
S4 = ElementTopology.QUAD4
MITC4 = ElementTopology.MITC4
S6 = ElementTopology.TRI6
S8 = ElementTopology.QUAD8
QSPT = ElementTopology.QSPT
