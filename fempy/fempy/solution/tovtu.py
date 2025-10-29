#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tovtk.py

Convert FEM geometry + solution to a single **VTK-HDF** file for ParaView
that contains the entire **time series** (all steps) in ONE container.

- If fields end with a numeric suffix (e.g. U0, U_1, Stress-12), we pack
  all steps into <stem>.vtkhdf with appropriate time values (float(idx)).
- If no time series (or --no_time), we write a single step with t=0.0.
- Output is placed in a new directory named after the geometry stem
  (or the --output stem), same behavior as before.

This version also:
- merges vector components (BASE_X/Y/Z and *_u_x/y/z) into 3-comp vectors,
- derives von Mises for 6-comp *stress* fields if missing,
- ensures 3-comp + magnitude for 3+ comp point fields (e.g., displacement/modes).

Requirements:
- A recent VTK build with VTKHDF writer support (vtkHDFWriter).
  If not found, the script raises a clear error message.

Author: Finn Eggers
Date  : 21.10.2025
"""

import argparse
import time
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

# local imports
from ..geometry import Geometry
from .solution import Solution


# ------------------------ logging helpers ------------------------

def log_info(msg: str) -> None:
    print(f"[INFO] {msg}")

def log_timing(t0: float, label: str) -> None:
    dt_ms = (time.time() - t0) * 1000.0
    print(f"[INFO] Finished {label:<60} [{dt_ms:8.2f} ms]")


# ------------------------ name/time utilities ------------------------

# Match anything that ends with a number (optionally preceded by _ or -)
_NUM_SUFFIX = re.compile(r"^(?P<base>.*?)(?:[_\-]?)(?P<idx>\d+)$", re.IGNORECASE)

def split_name(field_name: str):
    """
    Split a field name into (base, index) if it ends with a number.

    Examples:
        "U_T17"    -> ("U_T", 17)
        "U17"      -> ("U", 17)
        "TEMP_003" -> ("TEMP", 3)
        "PRESSURE" -> ("PRESSURE", None)
    """
    m = _NUM_SUFFIX.match(field_name)
    if not m:
        return field_name, None
    try:
        return m.group("base"), int(m.group("idx"))
    except ValueError:
        return field_name, None


def reorder_stress_if_needed(arr: np.ndarray) -> np.ndarray:
    """
    FEMaster: [XX,YY,ZZ, YZ,XZ,XY]  ->  VTK: [XX,YY,ZZ, XY,YZ,XZ]
    Only when arr has shape (N,6).
    """
    if arr is None:
        return None
    if arr.ndim == 2 and arr.shape[1] == 6:
        arr = arr.copy()
        arr[:, [3, 4, 5]] = arr[:, [5, 3, 4]]
    return arr


def numpy_to_vtk_array(data: np.ndarray, name: str) -> vtk.vtkDataArray:
    """
    Robust conversion of 1D/2D numpy arrays to vtkDataArray.
    - Deep copy to detach from numpy memory.
    - Explicitly set number of components for 2D arrays.
    """
    arr = np.asarray(data)
    if arr.ndim == 1:
        vtk_arr = numpy_to_vtk(arr, deep=1)
        vtk_arr.SetName(name)
        return vtk_arr
    elif arr.ndim == 2:
        n, m = arr.shape
        vtk_arr = numpy_to_vtk(arr.reshape(-1), deep=1)
        vtk_arr.SetNumberOfComponents(m)
        vtk_arr.SetName(name)
        return vtk_arr
    else:
        vtk_arr = numpy_to_vtk(arr.reshape(-1), deep=1)
        vtk_arr.SetName(name)
        return vtk_arr


# ------------------------ vector component combiner ------------------------

# Matches "..._X" / "...-Y" / "..._z" at the end
_END_COMP = re.compile(r"^(?P<base>.*?)[_\-](?P<comp>[xyzXYZ])$", re.IGNORECASE)
# Matches "..._u_x" / "...-u-z" (case-insensitive)
_EMBED_U = re.compile(r"^(?P<prefix>.*?)[_\-]u[_\-](?P<comp>[xyzXYZ])$", re.IGNORECASE)

def _combine_vector_components(fields: dict, n_points: int) -> dict:
    """
    Find scalar point fields with X/Y/Z components and merge into one 3-comp vector.
    Supports:
      - trailing components:  BASE_X, BASE_Y, BASE_Z
      - embedded 'u' style:   PREFIX_u_x, PREFIX_u_y, PREFIX_u_z
    Only merges arrays of length == n_points (point-data); cell-data is left as-is.

    If only X and Y exist, Z is filled with zeros.
    Original component fields are removed after merging.
    """
    groups = {}   # base -> {'x': arr, 'y': arr, 'z': arr, '_names': set(original_keys)}
    for name, arr in fields.items():
        if arr is None:
            continue
        a = np.asarray(arr)
        # Only merge point-wise scalars (1D or (N,1))
        if not ((a.ndim == 1 and len(a) == n_points) or (a.ndim == 2 and a.shape[0] == n_points and a.shape[1] == 1)):
            continue

        m1 = _END_COMP.match(name)
        m2 = _EMBED_U.match(name)

        if m1:
            base = m1.group('base')
            comp = m1.group('comp').lower()
        elif m2:
            base = m2.group('prefix')
            comp = m2.group('comp').lower()
        else:
            continue

        if comp not in ('x', 'y', 'z'):
            continue

        g = groups.setdefault(base, {'x': None, 'y': None, 'z': None, '_names': set()})
        g[comp] = a.reshape(-1) if a.ndim > 1 else a
        g['_names'].add(name)

    merged = {}
    to_skip = set()
    for base, g in groups.items():
        if g['x'] is None or g['y'] is None:
            continue
        x = g['x']
        y = g['y']
        z = g['z'] if g['z'] is not None else np.zeros_like(x)
        vec = np.column_stack([x, y, z])
        merged_name = base  # keep base as vector field name
        merged[merged_name] = vec
        # also helpful aliases if base looks like displacement
        if merged_name.lower() in ("u", "disp", "displacement"):
            merged["Displacement"] = vec
        to_skip |= g['_names']

    out = {}
    out.update(merged)
    for k, v in fields.items():
        if k not in to_skip:
            out[k] = v
    return out


# ------------------------ derived fields helpers ------------------------

def _is_stress_name(name: str) -> bool:
    return "stress" in name.lower()

def _ensure_derived_fields(fields: dict, n_points: int, n_cells: int) -> dict:
    """
    - For any 6-comp '*stress*' array: add '<base>_mises' if missing.
    - For any point field with >=3 comps: add '<base>_xyz' (first 3 comps) and '<base>_mag'.
    Avoids overwriting if names already present from Solution.
    """
    out = dict(fields)

    # pass 1: add von Mises for stress-like 6-comp arrays (point or cell)
    for name, arr in list(fields.items()):
        if arr is None:
            continue
        A = np.asarray(arr)
        if A.ndim == 2 and A.shape[1] >= 6 and _is_stress_name(name):
            base = name
            mises_name = f"{base}_mises"
            if mises_name not in out:
                sx, sy, sz, tyz, tzx, txy = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5]
                mises = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2
                                       + 6.0 * (txy ** 2 + tyz ** 2 + tzx ** 2)))
                out[mises_name] = mises

    # pass 2: for point fields with >=3 comps, ensure xyz & magnitude exist
    for name, arr in list(fields.items()):
        if arr is None:
            continue
        A = np.asarray(arr)
        if A.ndim == 2 and A.shape[0] == n_points and A.shape[1] >= 3:
            base = name
            xyz_name = f"{base}_xyz"
            mag_name = f"{base}_mag"
            if xyz_name not in out:
                out[xyz_name] = A[:, :3]
            if mag_name not in out:
                out[mag_name] = np.linalg.norm(A[:, :3], axis=1)

            # if it looks like displacement/mode, provide a canonical "Displacement"
            if base.lower() in ("u", "disp", "displacement") or "mode_shape" in base.lower() or "buckling_mode" in base.lower():
                out["Displacement"] = A[:, :3]

    return out


# ------------------------ VTK helpers ------------------------

VTK_CELL_CLASS = {
    'C3D4': vtk.vtkTetra,
    'C3D6': vtk.vtkWedge,
    'C3D8': vtk.vtkHexahedron,
    'C3D10': vtk.vtkQuadraticTetra,
    'C3D15': vtk.vtkQuadraticWedge,
    'C3D20': vtk.vtkQuadraticHexahedron,
    'C3D20R': vtk.vtkQuadraticHexahedron,
    'C2D3': vtk.vtkTriangle,
    'C2D4': vtk.vtkQuad,
    'C2D6': vtk.vtkQuadraticTriangle,
    'C2D8': vtk.vtkQuadraticQuad,
    'S4': vtk.vtkQuad,
    'S8': vtk.vtkQuadraticQuad,
    'S3': vtk.vtkTriangle,
    'S6': vtk.vtkQuadraticTriangle,
    "B33": vtk.vtkLine,
}


# ------------------------ HDF writer (single-file time series) ------------------------

class HDFSeriesWriter:
    """
    Thin convenience wrapper around vtkHDFWriter to write multiple time steps
    into a single .vtkhdf file. Requires a VTK build that provides vtkHDFWriter.
    """
    def __init__(self, hdf_path: Path):
        self.hdf_path = Path(hdf_path)
        self._first = True

        # Feature detection
        if not hasattr(vtk, "vtkHDFWriter"):
            raise RuntimeError(
                "vtkHDFWriter not found in this VTK build. "
                "Please use a recent VTK/ParaView (VTKHDF time-series support) or build VTK with HDF enabled."
            )

        self._writer = vtk.vtkHDFWriter()
        self._writer.SetFileName(str(self.hdf_path))

        # Optional compression knob (may not exist in all builds)
        if hasattr(self._writer, "SetCompressionLevel"):
            try:
                self._writer.SetCompressionLevel(5)
            except Exception:
                pass

    def write_step(self, grid: vtk.vtkUnstructuredGrid, time_value: float):
        """
        Append one step to the HDF file.
        - For static topology, we re-use the same connectivity each time.
        - Attaches time value when the API exposes it; otherwise the writer
          will still create the Steps group (recent VTK builds).
        """
        self._writer.SetInputData(grid)

        # Time metadata (API differs across VTK versions; try in order)
        if hasattr(self._writer, "SetTimeValue"):
            self._writer.SetTimeValue(float(time_value))
        elif hasattr(self._writer, "SetTimeStep"):
            # Fallback: step index; time values may be stored separately depending on build
            self._writer.SetTimeStep(int(time_value))

        # Append after first write (method name varies)
        if hasattr(self._writer, "SetAppend"):
            self._writer.SetAppend(not self._first)
        elif hasattr(self._writer, "SetWriteModeAppend"):
            self._writer.SetWriteModeAppend(not self._first)

        ok = self._writer.Write()
        if ok != 1:
            raise RuntimeError(f"vtkHDFWriter reported non-success when writing time={time_value}")

        self._first = False


# ------------------------ main converter class ------------------------

class Converter:
    def __init__(self, geometry_path, solution_path=None, output_filename=None,
                 force_no_time: bool = False):
        self.geometry = Geometry.read_input_deck(geometry_path)

        # solution
        if solution_path is None:
            solution_path = geometry_path.replace(".inp", ".res")

        log_info(f"Loading solution file '{solution_path}'")
        try:
            self.solution = Solution.open(solution_path)
            self.has_solution = True
            log_info(f"Solution file '{solution_path}' loaded successfully")
            print(self.solution.summary())
        except FileNotFoundError:
            log_info(f"Warning: Solution file '{solution_path}' not found. Geometry will be written without solution fields.")
            self.solution = None
            self.has_solution = False

        # Decide output directory & base stem (as before)
        if output_filename:
            out_path = Path(output_filename)
            self._stem = out_path.stem
            self._out_dir = (out_path.parent if str(out_path.parent) != "" else Path(".")) / self._stem
        else:
            inp_path = Path(geometry_path)
            parent = inp_path.parent if str(inp_path.parent) != "" else Path(".")
            self._stem = inp_path.stem
            self._out_dir = parent / self._stem

        self._out_dir.mkdir(parents=True, exist_ok=True)
        log_info(f"Output directory: {self._out_dir}")

        # base grid (geometry only); copied per frame when writing
        self._ugrid_geom = vtk.vtkUnstructuredGrid()

        # config
        self.force_no_time = bool(force_no_time)

    # ---------- public entry points ----------

    def convert(self):
        """
        Write a single .vtkhdf file containing either:
          - a single step (t=0.0), or
          - the full time series with one step per numeric suffix.
        """
        t0 = time.time()
        log_info("Begin conversion")
        self._generate_geometry()

        hdf_path = self._out_dir / f"{self._stem}.vtkhdf"

        # Geometry-only (no solution): single step at t=0
        if not self.has_solution:
            writer = HDFSeriesWriter(hdf_path)
            grid0 = vtk.vtkUnstructuredGrid()
            grid0.DeepCopy(self._ugrid_geom)
            writer.write_step(grid0, 0.0)
            log_timing(t0, "conversion (geometry-only, single .vtkhdf)")
            log_info(f"Wrote single-file time-series: {hdf_path}")
            return

        fields_dict = self.solution.list_fields_reduced(None)
        time_indices, grouped, static = self._group_fields_by_time(fields_dict)

        # Single snapshot mode (no suffixes or --no_time)
        if self.force_no_time or len(time_indices) == 0:
            log_info("No numeric-suffix fields detected or --no_time set → writing single step (t=0.0) to .vtkhdf")
            grid = self._grid_with_fields(static | self._flatten_time_groups(grouped))
            writer = HDFSeriesWriter(hdf_path)
            writer.write_step(grid, 0.0)
            log_timing(t0, "conversion (single snapshot, single .vtkhdf)")
            log_info(f"Wrote single-file time-series: {hdf_path}")
            return

        # Multi-time series → one .vtkhdf with all steps
        writer = HDFSeriesWriter(hdf_path)
        log_info(f"Begin writing single-file HDF time series: {hdf_path}")
        t_ts = time.time()

        # Use the numeric suffix as the time value (float(idx)), as before
        for k, idx in enumerate(time_indices):
            frame_fields = {}
            frame_fields.update(static)
            frame_fields.update(grouped[idx])

            grid_k = self._grid_with_fields(frame_fields)

            tval = float(idx)
            writer.write_step(grid_k, tval)
            log_info(f"  appended step {k+1}/{len(time_indices)}  idx={idx}  t={tval}")

        log_timing(t_ts, "writing HDF time series (single file)")
        log_timing(t0,  "conversion (HDF series)")
        log_info(f"Wrote single-file time-series: {hdf_path}")

    # ---------- geometry ----------

    def _generate_geometry(self):
        log_info("Begin geometry generation")
        t0 = time.time()
        self._set_points(self._ugrid_geom)
        self._set_cells(self._ugrid_geom)
        log_timing(t0, "geometry generation")

    def _set_points(self, ugrid: vtk.vtkUnstructuredGrid):
        log_info("Begin setting points (nodes)")
        t0 = time.time()
        points = vtk.vtkPoints()
        for node in self.geometry.nodes:
            points.InsertNextPoint(node if node is not None else (0.0, 0.0, 0.0))
        ugrid.SetPoints(points)
        log_timing(t0, "setting points (nodes)")

    def _set_cells(self, ugrid: vtk.vtkUnstructuredGrid):
        log_info("Begin setting cells (elements)")
        t0 = time.time()
        for element in self.geometry.elements:
            if element is None:
                ugrid.InsertNextCell(0, vtk.vtkIdList())
                continue
            cell_cls = VTK_CELL_CLASS.get(element.elem_type)
            if cell_cls is None:
                print(f"[WARNING] Unknown cell type: {element.elem_type}")
                ugrid.InsertNextCell(0, vtk.vtkIdList())
                continue
            cell = cell_cls()
            for j, nid in enumerate(element.node_ids):
                cell.GetPointIds().SetId(j, int(nid))
            ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        log_timing(t0, "setting cells (elements)")

    # ---------- solution handling ----------

    def _group_fields_by_time(self, fields_dict):
        grouped = defaultdict(dict)
        static = {}
        for name, func in fields_dict.items():
            base, idx = split_name(name)
            arr = func()
            if arr is None:
                continue
            if idx is None:
                static[base] = arr
            else:
                grouped[idx][base] = arr
        time_indices = sorted(grouped.keys())
        log_info(f"Detected {len(time_indices)} time steps from trailing numeric suffixes")
        return time_indices, grouped, static

    def _flatten_time_groups(self, grouped):
        flat = {}
        for idx, d in grouped.items():
            for base, arr in d.items():
                flat[f"{base}_{idx}"] = arr
        return flat

    def _grid_with_fields(self, fields_map):
        """
        Attach provided fields (numpy arrays) to either PointData or CellData.
        Also:
          - merges vector components to 3D vectors,
          - derives stresses' von Mises,
          - ensures 3-comp and magnitude for point vectors.
        """
        grid = vtk.vtkUnstructuredGrid()
        grid.DeepCopy(self._ugrid_geom)

        n_points = len(self.geometry.nodes)
        n_cells  = len(self.geometry.elements)

        # 1) merge component triplets to vectors
        fields_map = _combine_vector_components(fields_map, n_points)
        # 2) derive handy fields (mises, xyz, mag)
        fields_map = _ensure_derived_fields(fields_map, n_points, n_cells)

        for field_name, data in fields_map.items():
            if data is None:
                continue

            data = np.asarray(data)
            # Stress reordering just before writing
            if _is_stress_name(field_name) and data.ndim == 2 and data.shape[1] == 6:
                data = reorder_stress_if_needed(data)

            vtk_arr = numpy_to_vtk_array(data, field_name)

            if data.ndim == 1:
                # scalar point or cell?
                if len(data) == n_points:
                    grid.GetPointData().AddArray(vtk_arr)
                elif len(data) == n_cells:
                    grid.GetCellData().AddArray(vtk_arr)
                else:
                    print(f"[WARNING] Field '{field_name}' size {len(data)} does not match nodes ({n_points}) or elements ({n_cells}). Skipped.")
            elif data.ndim == 2:
                # vector/tensor per point or per cell?
                if data.shape[0] == n_points:
                    grid.GetPointData().AddArray(vtk_arr)
                elif data.shape[0] == n_cells:
                    grid.GetCellData().AddArray(vtk_arr)
                else:
                    print(f"[WARNING] Field '{field_name}' shape {data.shape} does not match nodes (*,{n_points}) or elements (*,{n_cells}). Skipped.")
            else:
                print(f"[WARNING] Field '{field_name}' has unsupported ndim={data.ndim}. Skipped.")

        return grid


# ------------------------ CLI ------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Convert geometry + solution to a single VTK-HDF file containing all time steps."
    )
    ap.add_argument("geometry_path", help="Path to the geometry input file (.inp).")
    ap.add_argument("--solution_path", default=None, help="Optional path to the solution file (.res).")
    ap.add_argument("--output", default=None, help="Base name for outputs. Directory will be <stem>/ next to this path; file is <stem>.vtkhdf.")
    ap.add_argument("--no_time", action="store_true", help="Force writing a single snapshot (no time series).")
    args = ap.parse_args()

    conv = Converter(
        geometry_path=args.geometry_path,
        solution_path=args.solution_path,
        output_filename=args.output,
        force_no_time=args.no_time,
    )
    conv.convert()


if __name__ == "__main__":
    main()
