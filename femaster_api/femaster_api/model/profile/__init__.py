"""Shared beam-profile definitions.

The package contains the reusable geometric beam-profile object and its named
project-level repository.  Profiles are deliberately separate from beam
sections: profiles describe section geometry/inertial properties, while
``BeamSection`` performs the material and element-region assignment.
"""

from .profile import Profile
from .profile_repository import ProfileRepository

__all__ = ["Profile", "ProfileRepository"]
