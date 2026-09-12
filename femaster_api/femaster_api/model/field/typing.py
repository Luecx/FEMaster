"""Address types used by model and result ``Field`` storage.

Editable model fields keep real FEM objects in every entity-address position.
``NODE`` fields are keyed by ``Node`` objects, ``ELEMENT`` fields by ``Element``
objects, and element-local domains use ``(Element, local_index, ...)`` tuples.
The local indices are intrinsic zero-based positions and therefore remain
integers.

Result files are different: RES/FRD addresses are solver output identifiers and
may be bare integers or qualified strings such as ``"bolt.17"``.  Those
serialized result addresses are retained deliberately because a ``Result`` does
not own the editable model object graph.
"""

from __future__ import annotations

from ..element.element import Element
from ..node.node import Node


ModelFieldKey = (
    Node
    | Element
    | tuple[Element, int]
    | tuple[Element, int, int]
)

FieldIndex = int | str
ResultFieldKey = FieldIndex | tuple[FieldIndex, ...]
FieldKey = ModelFieldKey | ResultFieldKey
