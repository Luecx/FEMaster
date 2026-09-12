"""Collector-owned FEMaster load definitions.

Every concrete load lives in its own ``load_*`` module and owns its native export
syntax.  ``LoadCollector`` provides the semantic grouping activated by analysis
steps, while ``LoadCollectorRepository`` provides project-level name lookup.
No second global repository duplicates individual load objects.
"""

from .load import Load
from .load_collector import LoadCollector
from .load_collector_repository import LoadCollectorRepository
from .load_inertial import InertialLoad
from .load_nodal_force import NodalForce
from .load_pressure import PressureLoad
from .load_surface_traction import SurfaceTraction
from .load_thermal import ThermalLoad
from .load_volume import VolumeLoad

__all__ = [
    "InertialLoad",
    "Load",
    "LoadCollector",
    "LoadCollectorRepository",
    "NodalForce",
    "PressureLoad",
    "SurfaceTraction",
    "ThermalLoad",
    "VolumeLoad",
]
