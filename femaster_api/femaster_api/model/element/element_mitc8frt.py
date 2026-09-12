"""Concrete ``MITC8FRT`` element definition.

Eight-node finite-rotation MITC shell element.  The class contains only the topology information that
distinguishes this native element type: its FEMaster ``TYPE`` token and required
number of connectivity nodes.  Persistent IDs and connectivity validation are
implemented once by the shared ``Element`` base class.

Keeping every concrete element in its own module makes the public element set
easy to extend without growing a central monolithic mesh file.
"""

from __future__ import annotations

from .element import Element


class MITC8FRT(Element):
    """Eight-node finite-rotation MITC shell element."""

    type_name = "MITC8FRT"
    node_count = 8
