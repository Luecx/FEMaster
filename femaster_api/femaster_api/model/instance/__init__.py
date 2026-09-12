"""Rigid assembly instances referencing reusable parts by semantic name.

Instances add only assembly placement; they never duplicate part topology.
``InstanceRepository`` owns deterministic instance order while ``Project`` owns
the surrounding assembly scope together with assembly regions and surfaces.
"""

from .instance import Instance
from .instance_repository import InstanceRepository

__all__ = ["Instance", "InstanceRepository"]
