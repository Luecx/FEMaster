"""Validation diagnostic objects."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Severity(Enum):
    ERROR = "error"
    WARNING = "warning"


@dataclass(frozen=True, slots=True)
class Diagnostic:
    severity: Severity
    message: str
    location: str | None = None


class ValidationError(Exception):
    """Raised when model validation finds at least one error."""


@dataclass
class Diagnostics:
    """Container for validation errors and warnings."""

    items: list[Diagnostic] = field(default_factory=list)

    def error(self, message: str, *, location: str | None = None) -> None:
        self.items.append(Diagnostic(Severity.ERROR, message, location))

    def warning(self, message: str, *, location: str | None = None) -> None:
        self.items.append(Diagnostic(Severity.WARNING, message, location))

    @property
    def errors(self) -> tuple[Diagnostic, ...]:
        return tuple(item for item in self.items if item.severity == Severity.ERROR)

    @property
    def warnings(self) -> tuple[Diagnostic, ...]:
        return tuple(item for item in self.items if item.severity == Severity.WARNING)

    @property
    def ok(self) -> bool:
        return not self.errors

    def raise_for_errors(self) -> None:
        if not self.errors:
            return
        lines = [f"{item.location + ': ' if item.location else ''}{item.message}" for item in self.errors]
        raise ValidationError("\n".join(lines))
