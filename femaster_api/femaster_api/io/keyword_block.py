"""Syntax-level keyword block used by the INP importer."""

from dataclasses import dataclass, field


@dataclass(slots=True)
class KeywordBlock:
    """One keyword line together with its following data rows."""

    name: str
    keys: dict[str, str] = field(default_factory=dict)
    data: list[list[str]] = field(default_factory=list)
    line: int = 0

    def key(self, name: str, default: str | None = None) -> str | None:
        return self.keys.get(name.upper().replace(" ", ""), default)
