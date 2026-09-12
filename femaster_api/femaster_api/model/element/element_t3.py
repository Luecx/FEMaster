"""Concrete ``T3`` element definition.

Native two-node truss element using the short T3 type name.  The class contains only the topology information that
distinguishes this native element type: its FEMaster ``TYPE`` token and required
number of connectivity nodes.  Persistent IDs and connectivity validation are
implemented once by the shared ``Element`` base class.

Keeping every concrete element in its own module makes the public element set
easy to extend without growing a central monolithic mesh file.
"""

from __future__ import annotations

from .element import Element


class T3(Element):
    """Native two-node truss element using the short T3 type name."""

    type_name = "T3"
    node_count = 2
