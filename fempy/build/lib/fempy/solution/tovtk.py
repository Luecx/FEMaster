"""
tovtk.py

Script for converting FEM geometry and solution data to VTK format for visualization in tools like ParaView.
Supports a variety of 2D and 3D element types and handles node-based and element-based solution fields.

Author: Finn Eggers
Date: 07.10.2024

Usage:
    python tovtk.py geometry.inp [--solution_path solution.res] [--output output.vtk]

Arguments:
    geometry_path: Path to the geometry input file (.inp).
    --solution_path: Optional path to the solution file (.res). If not provided, deduced from the input file.
    --output: Optional output VTK file name. If not provided, deduced from the input file by replacing '.inp' with '.vtk'.
"""

# library imports
import argparse
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import time

# imports from the same directory
from ..geometry import Geometry
from .solution import Solution


def log_info(message):
    """Prints an info message with the [INFO] tag."""
    print(f"[INFO] {message}")

def log_timing(start_time, task_name):
    """Prints the time taken for a task."""
    elapsed_time = (time.time() - start_time) * 1000
    print(f"[INFO] Finished {task_name:<60} [{elapsed_time:8.2f} ms]")


class Converter:
    def __init__(self, geometry_path, solution_path=None, output_filename=None):
        self.geometry = Geometry.read_input_deck(geometry_path)

        # Deduce solution path if not provided
        if solution_path is None:
            solution_path = geometry_path.replace(".inp", ".res")

        # Try loading the solution file
        log_info(f"Loading solution file '{solution_path}'")
        try:
            self.solution = Solution.open(solution_path)
            self.has_solution = True
            log_info(f"Solution file '{solution_path}' loaded successfully")
        except FileNotFoundError:
            log_info(f"Warning: Solution file '{solution_path}' not found. Geometry will be written without solution fields.")
            self.solution = None
            self.has_solution = False

        # Deduce output filename if not provided
        self.output_filename = output_filename or geometry_path.replace(".inp", ".vtk")
        self._ugrid = vtk.vtkUnstructuredGrid()

    def convert(self):
        """Main function to convert geometry and optional solution data to a VTK file."""
        log_info("Begin conversion process")
        start_time = time.time()
        self._generate_geometry()
        if self.has_solution:
            self._add_solution_fields()
        self._write_vtk()
        log_timing(start_time, "conversion process")

    def _generate_geometry(self):
        """Generates the VTK unstructured grid from the geometry data."""
        log_info("Begin geometry generation")
        start_time = time.time()
        self._set_points()
        self._set_cells()
        log_timing(start_time, "geometry generation")

    def _set_points(self):
        """Set up the points (nodes) in the unstructured grid."""
        log_info("Begin setting points (nodes)")
        start_time = time.time()
        points = vtk.vtkPoints()
        for node in self.geometry.nodes:
            points.InsertNextPoint(node if node else [0, 0, 0])
        self._ugrid.SetPoints(points)
        log_timing(start_time, "setting points (nodes)")

    def _set_cells(self):
        """Define the cells (elements) in the unstructured grid."""
        log_info("Begin setting cells (elements)")
        start_time = time.time()
        for element in self.geometry.elements:
            cell_type, cell = self._create_vtk_cell(element)
            if cell:
                self._ugrid.InsertNextCell(cell_type, cell.GetPointIds())
            else:
                self._ugrid.InsertNextCell(0, vtk.vtkIdList())
        log_timing(start_time, "setting cells (elements)")

    def _create_vtk_cell(self, element):
        """Create a VTK cell based on the element type."""
        if element is None:
            return None, None
        cell_map = {
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

        cell_class = cell_map.get(element.elem_type)
        if cell_class is None:
            print (f'[WARNING] Unknown cell type: {element.elem_type}')
            return None, None
        cell = cell_class()
        for j, node_id in enumerate(element.node_ids):
            cell.GetPointIds().SetId(j, node_id)
        return cell.GetCellType(), cell

    def _add_solution_fields(self):
        """Add solution fields to the VTK file if available."""
        log_info("Begin adding solution fields")
        start_time = time.time()
        fields_dict = self.solution.list_fields_reduced(None)
        for field_name, field_func in fields_dict.items():
            field_data = field_func().copy()
            vtk_data = numpy_to_vtk(field_data)
            vtk_data.SetName(field_name)

            log_info(f"   Adding field '{field_name}' with shape {field_data.shape}")

            # if its 6 columns, reorder them to match VTK's ordering
            # FEMaster uses XX, YY, ZZ, YZ, XZ, XY
            # VTK      uses XX, YY, ZZ, XY, YZ, XZ
            if len(field_data.shape) > 1 and field_data.shape[1] == 6:
                field_data[:, [3, 4, 5]] = field_data[:, [5, 3, 4]]

            if len(field_data) == len(self.geometry.nodes):
                self._ugrid.GetPointData().AddArray(vtk_data)
            elif len(field_data) == len(self.geometry.elements):
                self._ugrid.GetCellData().AddArray(vtk_data)
            else:
                print(f"[WARNING] Field '{field_name}' size does not match nodes or elements.")
        log_timing(start_time, "adding solution fields")

    def _write_vtk(self):
        """Write the generated unstructured grid to a VTK file."""
        log_info(f"Begin writing VTK file: {self.output_filename}")
        start_time = time.time()
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(self.output_filename)
        writer.SetFileTypeToBinary()  # Set the file type to binary
        writer.SetInputData(self._ugrid)
        writer.Write()
        log_timing(start_time, "writing VTK file")


def main():
    parser = argparse.ArgumentParser(description="Convert geometry and solution data to VTK format.")
    parser.add_argument("geometry_path", help="Path to the geometry input file (.inp).")
    parser.add_argument("--solution_path", help="Optional path to the solution file (.res).", default=None)
    parser.add_argument("--output", help="Output VTK file name.", default=None)

    args = parser.parse_args()

    # Create the converter instance and perform the conversion
    converter = Converter(geometry_path=args.geometry_path, solution_path=args.solution_path, output_filename=args.output)
    converter.convert()


if __name__ == '__main__':
    main()
