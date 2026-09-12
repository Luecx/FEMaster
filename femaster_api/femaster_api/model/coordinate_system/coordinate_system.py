"""Base coordinate-system definition."""

from ..common.named_object import NamedObject


class CoordinateSystem(NamedObject):
    """Base named global coordinate system."""

    def export(self) -> str:
        raise NotImplementedError
