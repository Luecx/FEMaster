"""FEMaster input-deck reader."""

from __future__ import annotations

from pathlib import Path

from femaster_api.model.model import Model
from femaster_api.importers.femaster_deck import load_model_from_inp


class FEMasterReader:
    """Parse FEMaster `.inp` decks into the Python API model."""

    def read(self, path: str | Path) -> Model:
        return load_model_from_inp(path)
