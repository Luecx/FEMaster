"""Input and result readers for the FEMaster Python API."""

from __future__ import annotations

from pathlib import Path

from .frd import FrdReader
from .inp import InpReader, KeywordBlock
from .res import ResReader


def read_input(path: str | Path):
    """Read a FEMaster/Abaqus-like INP deck into a Project."""

    return InpReader().read(path)


def read_result(path: str | Path):
    """Read a supported FEMaster result file by filename extension."""

    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".res":
        return ResReader().read(path)
    if suffix == ".frd":
        return FrdReader().read(path)

    raise ValueError(f"unsupported result format: {path.suffix or '<none>'}")


__all__ = [
    "FrdReader",
    "InpReader",
    "KeywordBlock",
    "ResReader",
    "read_input",
    "read_result",
]
