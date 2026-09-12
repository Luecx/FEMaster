"""Input and result importers."""

from __future__ import annotations

from pathlib import Path

from .frd_importer import FrdImporter
from .inp_importer import InpImporter
from .keyword_block import KeywordBlock
from .res_importer import ResImporter


def import_input(path: str | Path):
    """Import a FEMaster/Abaqus-like input deck."""

    return InpImporter().import_file(path)


def import_result(path: str | Path):
    """Import a supported FEMaster result file by extension."""

    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".res":
        return ResImporter().import_file(path)
    if suffix == ".frd":
        return FrdImporter().import_file(path)

    raise ValueError(
        f"unsupported result format: {path.suffix or '<none>'}"
    )


__all__ = [
    "FrdImporter",
    "InpImporter",
    "KeywordBlock",
    "ResImporter",
    "import_input",
    "import_result",
]
