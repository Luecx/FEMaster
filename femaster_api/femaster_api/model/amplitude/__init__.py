"""Global amplitude curves referenced by time-dependent definitions.

An ``Amplitude`` is a named sequence of scalar time/value samples with an
explicit interpolation rule.  The project-level repository makes those curves
reusable across loads without embedding time histories in every load object.
"""

from .amplitude import Amplitude
from .amplitude_interpolation import AmplitudeInterpolation
from .amplitude_repository import AmplitudeRepository

__all__ = ["Amplitude", "AmplitudeInterpolation", "AmplitudeRepository"]
