"""Reusable support definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.orientations import Orientation

Target = str | int | object
SupportValues = tuple[float | None, float | None, float | None, float | None, float | None, float | None]


@dataclass(frozen=True, slots=True)
class Support:
    target: Target
    values: SupportValues
    orientation: Orientation | None = None
