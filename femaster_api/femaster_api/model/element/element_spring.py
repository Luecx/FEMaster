"""Concrete ``SPRING1`` element definition.

One-node SPRING1 topology for a ground spring.  The class contains only the topology information that
distinguishes this native element type: its FEMaster ``TYPE`` token and required
number of connectivity nodes.  Persistent IDs and connectivity validation are
implemented once by the shared ``Element`` base class.

Keeping every concrete element in its own module makes the public element set
easy to extend without growing a central monolithic mesh file.
"""

from __future__ import annotations

from .element import Element


class SpringElement(Element):
    """One-node SPRING1 topology for a ground spring."""

    type_name = "SPRING1"
    node_count = 1
