"""Concrete ``B33`` element definition.

Two-node beam element using the FEMaster B33 formulation.  The class contains only the topology information that
distinguishes this native element type: its FEMaster ``TYPE`` token and required
number of connectivity nodes.  Persistent IDs and connectivity validation are
implemented once by the shared ``Element`` base class.

Keeping every concrete element in its own module makes the public element set
easy to extend without growing a central monolithic mesh file.
"""

from __future__ import annotations

from .element import Element


class B33(Element):
    """Two-node beam element using the FEMaster B33 formulation."""

    type_name = "B33"
    node_count = 2
