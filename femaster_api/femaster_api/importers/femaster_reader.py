"""Minimal FEMaster input-deck reader scaffold.

The public package is primarily a modeling API plus exporter. Deck import is
kept deliberately small here because a robust importer should be implemented
against the FEMaster DSL registry rather than by duplicating parser behavior in
Python.
"""

from __future__ import annotations

from pathlib import Path

from femaster_api.model.model import Model


class FEMasterReader:
    """Placeholder for future full deck import support."""

    def read(self, path: str | Path) -> Model:
        raise NotImplementedError(
            "FEMasterReader is a scaffold. Use FEMasterWriter for export and "
            "ResultReader for solution extraction until deck import is implemented."
        )
