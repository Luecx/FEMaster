"""Load model objects."""

from .amplitude_interpolation import AmplitudeInterpolation
from .amplitude import Amplitude
from .inertial_load import InertialLoad
from .load_collector import LoadCollector
from .load_collector_repository import LoadCollectorRepository
from .load_repository import LoadRepository
from .nodal_force import NodalForce, NodalTarget
from .pressure_load import PressureLoad, SurfaceTarget
from .surface_traction import SurfaceTraction
from .thermal_load import ThermalLoad
from .volume_load import ElementTarget, VolumeLoad

LoadEntry = NodalForce | SurfaceTraction | PressureLoad | VolumeLoad | ThermalLoad | InertialLoad

__all__ = [
    "Amplitude",
    "AmplitudeInterpolation",
    "InertialLoad",
    "LoadCollector",
    "LoadCollectorRepository",
    "LoadEntry",
    "LoadRepository",
    "ElementTarget",
    "NodalForce",
    "NodalTarget",
    "PressureLoad",
    "SurfaceTarget",
    "SurfaceTraction",
    "ThermalLoad",
    "VolumeLoad",
]
