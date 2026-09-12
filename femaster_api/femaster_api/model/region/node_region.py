"""Typed NodeRegion."""

from .region import Region


class NodeRegion(Region):
    """FEMaster node region."""

    keyword_name = "NSET"
