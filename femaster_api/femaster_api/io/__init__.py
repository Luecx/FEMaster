"""Input and result importers for the FEMaster Python API."""

from __future__ import annotations

from pathlib import Path

from .frd import FrdReader as _FrdReader
from .inp import InpReader as _InpReader, KeywordBlock
from .res import ResReader as _ResReader


class InpImporter:
    """Import FEMaster/Abaqus-like INP decks into Project objects."""

    def __init__(self) -> None:
        self._impl = _InpReader()

    def import_file(self, path: str | Path):
        """Import one UTF-8 input file."""

        return self._impl.read(path)

    def import_text(self, text: str):
        """Import one input deck from text."""

        return self._impl.parse(text)


class ResImporter:
    """Import native FEMaster RES result files."""

    def __init__(self) -> None:
        self._impl = _ResReader()

    def import_file(self, path: str | Path):
        return self._impl.read(path)

    def import_text(self, text: str):
        return self._impl.parse(text)


class FrdImporter:
    """Import FEMaster-generated CalculiX/CGX FRD result files."""

    def __init__(self) -> None:
        self._impl = _FrdReader()

    def import_file(self, path: str | Path):
        return self._impl.read(path)

    def import_text(self, text: str):
        return self._impl.parse(text)


def import_input(path: str | Path):
    """Import a FEMaster/Abaqus-like INP deck into a Project."""

    return InpImporter().import_file(path)


def import_result(path: str | Path):
    """Import a supported FEMaster result file by filename extension."""

    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".res":
        return ResImporter().import_file(path)
    if suffix == ".frd":
        return FrdImporter().import_file(path)

    raise ValueError(f"unsupported result format: {path.suffix or '<none>'}")


__all__ = [
    "FrdImporter",
    "InpImporter",
    "KeywordBlock",
    "ResImporter",
    "import_input",
    "import_result",
]
