"""Base class for named FEMaster surface definitions.

Surfaces are geometric/topological definitions, not region membership lists.
Concrete subclasses describe how a surface is constructed from elements or
nodes.  Both part-local and assembly-level surfaces use the same Python classes;
their owner determines the scope in which export places the ``*SURFACE`` block.

The base class intentionally contains no serialization branch on concrete type.
Each supported surface construction owns its own ``export`` implementation.
"""

from __future__ import annotations

from ..common.named_object import NamedObject


class Surface(NamedObject):
    """Base class for one named surface definition."""

    def export(self) -> str:
        """Export this surface to its native ``*SURFACE`` representation."""

        raise NotImplementedError
