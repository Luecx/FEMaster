"""Central object container for FEMaster models."""

from __future__ import annotations

from .analysis import LoadcaseRepository
from .constraints import ConstraintRepository
from .elements import ElementRepository
from .features import FeatureRepository
from .fields import FieldRepository
from .loads import LoadCollectorRepository, LoadRepository
from .materials import MaterialRepository
from .nodes import NodeRepository
from .orientations import OrientationRepository
from .sections import SectionRepository
from .sets import ElementSetRepository, NodeSetRepository, SurfaceSetRepository
from .supports import SupportCollectorRepository, SupportRepository
from .surfaces import SurfaceRepository


class Model:
    """Neutral finite-element model container.

    The model is intentionally a collection of object repositories. It does
    not render FEMaster text and it does not assign FEMaster IDs.
    """

    def __init__(self, name: str = "model") -> None:
        self.name               = name
        self.nodes              = NodeRepository()
        self.elements           = ElementRepository()
        self.node_sets          = NodeSetRepository()
        self.element_sets       = ElementSetRepository()
        self.surface_sets       = SurfaceSetRepository()
        self.surfaces           = SurfaceRepository()
        self.orientations       = OrientationRepository()
        self.materials          = MaterialRepository()
        self.sections           = SectionRepository()
        self.fields             = FieldRepository()
        self.features           = FeatureRepository()
        self.loads              = LoadRepository()
        self.load_collectors    = LoadCollectorRepository()
        self.supports = SupportRepository()
        self.support_collectors = SupportCollectorRepository()
        self.constraints        = ConstraintRepository()
        self.loadcases          = LoadcaseRepository()

    def validate(self, *, raise_on_error: bool = False):
        """Validate object ownership and cross-repository references."""

        from femaster_api.validation.checks import validate_model

        diagnostics = validate_model(self)
        if raise_on_error:
            diagnostics.raise_for_errors()
        return diagnostics
