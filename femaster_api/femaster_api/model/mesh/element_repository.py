"""Repository for part-local finite elements."""

from ..common.format import block, keyword
from ..common.id_repository import IdRepository
from .element import Element


class ElementRepository(IdRepository[Element]):
    """Element repository keyed by semantic FEMaster element ID."""

    def export(self) -> str:
        groups: dict[str, list[Element]] = {}
        for element in self:
            groups.setdefault(element.type_name, []).append(element)

        result: list[str] = []
        for type_name, elements in groups.items():
            result.append(
                block([
                    keyword("ELEMENT", TYPE=type_name),
                    *(element.export() for element in elements),
                ])
            )
        return "\n\n".join(result)
