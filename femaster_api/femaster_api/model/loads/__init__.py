"""Load model objects."""

from .amplitude_interpolation import AmplitudeInterpolation
from .amplitude import Amplitude
from .inertial_load import InertialLoad
from .load_collector import LoadCollector
from .load_collector_repository import LoadCollectorRepository
from .load_repository import LoadRepository
from .nodal_force import NodalForce
from .pressure_load import PressureLoad
from .surface_traction import SurfaceTraction
from .thermal_load import ThermalLoad
from .volume_load import VolumeLoad

LoadEntry = NodalForce | SurfaceTraction | PressureLoad | VolumeLoad | ThermalLoad | InertialLoad

__all__ = [
    "Amplitude",
    "AmplitudeInterpolation",
    "InertialLoad",
    "LoadCollector",
    "LoadCollectorRepository",
    "LoadEntry",
    "LoadRepository",
    "NodalForce",
    "PressureLoad",
    "SurfaceTraction",
    "ThermalLoad",
    "VolumeLoad",
]
