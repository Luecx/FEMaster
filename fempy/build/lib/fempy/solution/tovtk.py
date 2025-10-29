# file: tovtk.py
"""
tovtk.py  —  Batch converter (multi-.inp → .vtk in same directories)

Converts one or more FEM .inp + optional .res files to VTK (legacy .vtk).
Each output .vtk is written alongside its corresponding .inp.

Author: Finn Eggers
Date: 29.10.2025

Examples:
  python tovtk.py mesh1.inp mesh2.inp
  python tovtk.py ./cases --recursive
  python tovtk.py ./cases --solution_ext .dat --workers 8
"""

from __future__ import annotations
import argparse, time, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

# local imports
from ..geometry import Geometry
from .solution import Solution


# -----------------------------------------------------------------------------#
def log_info(msg): print(f"[INFO] {msg}")
def log_warn(msg): print(f"[WARNING] {msg}")
def log_timing(start, name): print(f"[INFO] Finished {name:<55} [{(time.time()-start)*1000:8.2f} ms]")


# -----------------------------------------------------------------------------#
class Converter:
    def __init__(self, geometry_path: Path, solution_path: Path | None, output_filename: Path):
        self.geometry_path = geometry_path
        self.output_filename = output_filename
        self.geometry = Geometry.read_input_deck(str(geometry_path))

        # Try to load optional solution
        self.solution = None
        self.has_solution = False
        if solution_path is not None:
            log_info(f"Loading solution '{solution_path}'")
            try:
                self.solution = Solution.open(str(solution_path))
                self.has_solution = True
            except FileNotFoundError:
                log_warn(f"Solution file '{solution_path}' not found – geometry only.")

        self._ugrid = vtk.vtkUnstructuredGrid()

    # -------------------------------------------------------------------------
    def convert(self):
        log_info(f"Converting '{self.geometry_path.name}'")
        t0 = time.time()
        self._generate_geometry()
        if self.has_solution:
            self._add_solution_fields()
        self._write_vtk()
        log_timing(t0, f"conversion '{self.geometry_path.name}'")

    def _generate_geometry(self):
        t0 = time.time()
        self._set_points()
        self._set_cells()
        log_timing(t0, "geometry generation")

    def _set_points(self):
        pts = vtk.vtkPoints()
        for n in self.geometry.nodes:
            pts.InsertNextPoint(n if n is not None else (0,0,0))
        self._ugrid.SetPoints(pts)

    def _set_cells(self):
        cell_map = {
            'C3D4': vtk.vtkTetra, 'C3D6': vtk.vtkWedge, 'C3D8': vtk.vtkHexahedron,
            'C3D10': vtk.vtkQuadraticTetra, 'C3D15': vtk.vtkQuadraticWedge,
            'C3D20': vtk.vtkQuadraticHexahedron, 'C3D20R': vtk.vtkQuadraticHexahedron,
            'C2D3': vtk.vtkTriangle, 'C2D4': vtk.vtkQuad, 'C2D6': vtk.vtkQuadraticTriangle,
            'C2D8': vtk.vtkQuadraticQuad, 'S3': vtk.vtkTriangle, 'S4': vtk.vtkQuad,
            'S6': vtk.vtkQuadraticTriangle, 'S8': vtk.vtkQuadraticQuad, 'B33': vtk.vtkLine,
        }
        for e in self.geometry.elements:
            if e is None:
                self._ugrid.InsertNextCell(0, vtk.vtkIdList()); continue
            cls = cell_map.get(e.elem_type)
            if not cls: log_warn(f"Unknown elem {e.elem_type}"); continue
            cell = cls()
            for j,nid in enumerate(e.node_ids): cell.GetPointIds().SetId(j,int(nid))
            self._ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

    def _add_solution_fields(self):
        fields = self.solution.list_fields_reduced(None)
        n_nodes, n_elems = len(self.geometry.nodes), len(self.geometry.elements)
        for name, func in fields.items():
            data = np.asarray(func())
            if data.ndim==2 and data.shape[1]==6:
                data = data.copy(); data[:,[3,4,5]] = data[:,[5,3,4]]  # reorder
            vtk_data = numpy_to_vtk(np.ascontiguousarray(data), deep=True)
            if data.ndim==2: vtk_data.SetNumberOfComponents(data.shape[1])
            vtk_data.SetName(name)
            if len(data)==n_nodes: self._ugrid.GetPointData().AddArray(vtk_data)
            elif len(data)==n_elems: self._ugrid.GetCellData().AddArray(vtk_data)
            else: log_warn(f"'{name}' size {len(data)} mismatch")

    def _write_vtk(self):
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(str(self.output_filename))
        writer.SetFileTypeToBinary()
        writer.SetInputData(self._ugrid)
        writer.Write()
        log_info(f"→ {self.output_filename}")


# -----------------------------------------------------------------------------#
def gather_inp_files(inputs, recursive=False):
    files=[]
    for p in inputs:
        pth=Path(p)
        if pth.is_dir():
            files+=sorted(pth.rglob("*.inp") if recursive else pth.glob("*.inp"))
        elif pth.is_file() and pth.suffix.lower()==".inp":
            files.append(pth)
    seen=set(); uniq=[]
    for f in files:
        if f not in seen: seen.add(f); uniq.append(f)
    return uniq

def deduce_paths(geometry_path:Path, solution_ext:str)->tuple[Path,Path]:
    stem=geometry_path.with_suffix("")
    return geometry_path.with_suffix(solution_ext), geometry_path.with_suffix(".vtk")

def convert_one(geometry_path:Path, solution_ext:str)->tuple[bool,str]:
    try:
        sol_path,out_path=deduce_paths(geometry_path,solution_ext)
        conv=Converter(geometry_path,sol_path,out_path)
        conv.convert()
        return True,str(out_path)
    except Exception as e:
        return False,f"{geometry_path}: {e}"


# -----------------------------------------------------------------------------#
def main():
    ap=argparse.ArgumentParser(description="Convert one or many FEM .inp files to .vtk (written next to input).")
    ap.add_argument("inputs",nargs="+",help="One or more .inp files and/or directories.")
    ap.add_argument("--solution_ext",type=str,default=".res",help="Extension of solution file (default: .res).")
    ap.add_argument("--recursive",action="store_true",help="Recurse into directories.")
    ap.add_argument("--workers",type=int,default=1,help="Parallel workers.")
    args=ap.parse_args()

    inp_files=gather_inp_files(args.inputs,args.recursive)
    if not inp_files:
        log_warn("No .inp files found."); return

    log_info(f"Found {len(inp_files)} file(s). Output beside each .inp.")
    ok=fail=0

    if args.workers<=1:
        for f in inp_files:
            success,msg=convert_one(f,args.solution_ext)
            ok+=success; fail+=not success
            log_info(f"OK  -> {msg}" if success else f"ERR -> {msg}")
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs=[ex.submit(convert_one,f,args.solution_ext) for f in inp_files]
            for fut in as_completed(futs):
                success,msg=fut.result()
                ok+=success; fail+=not success
                log_info(f"OK  -> {msg}" if success else f"ERR -> {msg}")

    print(f"\nDone. ok={ok}, failed={fail}.")


if __name__=="__main__":
    main()
