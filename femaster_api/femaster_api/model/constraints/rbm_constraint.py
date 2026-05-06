"""Rigid-body-motion suppression constraint."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import EntitySet


@dataclass(frozen=True, slots=True)
class RBMConstraint:
    element_set: EntitySet
