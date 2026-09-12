"""Typed ElementRegion."""

from .region import Region


class ElementRegion(Region):
    """FEMaster element region."""

    keyword_name = "ELSET"
