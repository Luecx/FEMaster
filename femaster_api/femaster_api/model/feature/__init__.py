"""Assembly-level non-topological model features.

Features are solver/model definitions that act on compiled entities without
being mesh topology, sections, loads or kinematic constraints.  Each concrete
feature remains isolated in its own module and is owned in deterministic order
by ``FeatureRepository``.
"""

from .feature import Feature
from .feature_point_mass import PointMass
from .feature_repository import FeatureRepository

__all__ = ["Feature", "FeatureRepository", "PointMass"]
