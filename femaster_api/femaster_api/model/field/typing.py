"""Semantic key aliases used by sparse ``Field`` storage.

A node/element field usually uses one integer or qualified string key.
Element-location fields use tuples such as ``("bolt.17", local_ip)`` and
material-point fields use ``("bolt.17", local_ip, local_mp)``.
"""

from __future__ import annotations

FieldIndex = int | str
FieldKey = FieldIndex | tuple[FieldIndex, ...]
