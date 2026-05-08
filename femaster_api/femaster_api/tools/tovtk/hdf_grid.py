"""VTKHDF unstructured-grid writer for FEMaster models."""

from __future__ import annotations

from pathlib import Path

from femaster_api.model import Model

from .mesh_export import VtkMesh, build_vtk_mesh


class VtkHdfWriter:
    """Write a model mesh to a `.vtkhdf` unstructured-grid file."""

    def __init__(self, model: Model) -> None:
        self.model = model

    def write(self, path: str | Path) -> None:
        """Write the model mesh using the VTKHDF layout."""

        mesh = build_vtk_mesh(self.model)
        _write_hdf(mesh, path)


def write_vtkhdf(model: Model, path: str | Path) -> None:
    """Convenience function for writing a model mesh to `.vtkhdf`."""

    VtkHdfWriter(model).write(path)


def _write_hdf(mesh: VtkMesh, path: str | Path) -> None:
    try:
        import h5py
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("VTKHDF export requires optional dependencies: h5py and numpy") from exc

    points = np.asarray(mesh.points, dtype=np.float64)
    connectivity = np.asarray(mesh.connectivity, dtype=np.int64)
    offsets = np.asarray((0, *mesh.offsets), dtype=np.int64)
    cell_types = np.asarray(mesh.cell_types, dtype=np.uint8)

    with h5py.File(path, "w") as h5:
        root = h5.create_group("VTKHDF")
        root.attrs["Version"] = (2, 2)
        root.attrs["Type"] = "UnstructuredGrid"
        root.create_dataset("Points", data=points)
        root.create_dataset("Connectivity", data=connectivity)
        root.create_dataset("Offsets", data=offsets)
        root.create_dataset("Types", data=cell_types)
