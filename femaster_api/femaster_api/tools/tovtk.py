"""Compatibility wrapper for VTK export."""

from __future__ import annotations

from femaster_api.tools.vtk_export import tovtk, write_vtk, write_vtkhdf
from femaster_api.tools.vtk_export.cli import main

__all__ = ["main", "tovtk", "write_vtk", "write_vtkhdf"]


if __name__ == "__main__":
    main()
