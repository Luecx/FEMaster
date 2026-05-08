"""Rigid-body-motion suppression constraint."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import ElementSet


@dataclass(frozen=True, slots=True)
class RBMConstraint:
    """Rigid-body-motion suppression over an element set."""

    element_set: ElementSet
