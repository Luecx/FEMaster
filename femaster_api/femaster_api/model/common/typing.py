"""Shared lightweight type aliases used by the FEMaster model.

The aliases in this module describe semantic references without introducing
wrapper classes solely for typing.  Integer entity references denote local
FEMaster IDs, while strings may denote named regions or instance-qualified
identifiers such as ``"bolt.17"``.
"""

from __future__ import annotations

EntityReference = int | str
NodeReference = int | str
ElementReference = int | str
