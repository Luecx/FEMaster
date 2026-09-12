"""Repository of global coordinate systems."""

from ..common.named_repository import NamedRepository
from .coordinate_system import CoordinateSystem


class CoordinateSystemRepository(NamedRepository[CoordinateSystem]):
    """Named global coordinate-system repository."""

    def export(self) -> str:
        return "\n\n".join(item.export() for item in self)
