"""Base material-backed section."""

from .section import Section


class MaterialSection(Section):
    """Base section referencing a global material."""

    def __init__(self, name: str, element_region: str, material: str) -> None:
        super().__init__(name, element_region)
        self.material = str(material)
