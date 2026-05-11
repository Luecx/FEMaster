"""Central model container for the independent FEMaster Python API."""

from __future__ import annotations

from .analysis import StepRepository
from .constraints import ConstraintRepository
from .elements import ElementRepository
from .features import FeatureRepository
from .fields import FieldRepository
from .loads import LoadCollectorRepository, LoadRepository
from .materials import MaterialRepository
from .nodes import NodeRepository
from .orientations import OrientationRepository
from .profiles import ProfileRepository
from .sections import SectionRepository
from .sets import SetRepository
from .supports import SupportCollectorRepository, SupportRepository
from .surfaces import SurfaceRepository


class Model:
    """Neutral finite-element model container.

    The model does not render FEMaster text itself. Use
    `femaster_api.export.FEMasterWriter` for deck generation.
    """

    def __init__(self, name: str = "model") -> None:
        self.name = name
        self.nodes = NodeRepository()
        self.elements = ElementRepository()
        self.sets = SetRepository()
        self.surfaces = SurfaceRepository()
        self.orientations = OrientationRepository()
        self.profiles = ProfileRepository()
        self.materials = MaterialRepository()
        self.sections = SectionRepository()
        self.fields = FieldRepository()
        self.features = FeatureRepository()
        self.loads = LoadRepository()
        self.load_collectors = LoadCollectorRepository()
        self.supports = SupportRepository()
        self.support_collectors = SupportCollectorRepository()
        self.constraints = ConstraintRepository()
        self.steps = StepRepository()

    def validate(self, *, raise_on_error: bool = False):
        from femaster_api.validation.checks import validate_model

        diagnostics = validate_model(self)
        if raise_on_error:
            diagnostics.raise_for_errors()
        return diagnostics
