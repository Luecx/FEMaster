"""Concrete ``ROTARYI`` element definition.

One-node ROTARYI topology for concentrated rotational inertia.  The class contains only the topology information that
distinguishes this native element type: its FEMaster ``TYPE`` token and required
number of connectivity nodes.  Persistent IDs and connectivity validation are
implemented once by the shared ``Element`` base class.

Keeping every concrete element in its own module makes the public element set
easy to extend without growing a central monolithic mesh file.
"""

from __future__ import annotations

from .element import Element


class RotaryInertiaElement(Element):
    """One-node ROTARYI topology for concentrated rotational inertia."""

    type_name = "ROTARYI"
    node_count = 1
