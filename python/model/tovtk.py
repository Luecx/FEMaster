import vtk
import numpy as np
import argparse
from vtk.util.numpy_support import numpy_to_vtk
from .geometry import Geometry
from .solution import Solution

class Convert:
    def __init__(self, geometry_path, solution_path=None, output_filename=None):
        self.geometry = Geometry.read_input_deck(geometry_path)

        # Deduce solution path if not provided
        if solution_path is None:
            solution_path = geometry_path.replace(".inp", ".res")
        self.solution = Solution.open(solution_path)

        # Deduce output filename if not provided
        if output_filename is None:
            output_filename = geometry_path.replace(".inp", ".vtk")

        self.output_filename = output_filename
        self.ugrid = None

    def _generate_ugrid(self):
        # Create the unstructured grid object
        self.ugrid = vtk.vtkUnstructuredGrid()

        # Set the points (geometry nodes)
        points = vtk.vtkPoints()
        for node in self.geometry.nodes:
            if node:
                points.InsertNextPoint(node)
            else:
                points.InsertNextPoint([0, 0, 0])
        self.ugrid.SetPoints(points)

        print(self.geometry)

        # Set the cells (geometry elements)
        for element in self.geometry.elements:
            if element:
                cellType = element.elem_type
                nodes = element.node_ids
                if cellType == 'C3D4':
                    cell = vtk.vtkTetra()
                elif cellType == 'C3D6':
                    cell = vtk.vtkWedge()
                elif cellType == 'C3D8':
                    cell = vtk.vtkHexahedron()
                elif cellType == 'C3D10':
                    cell = vtk.vtkQuadraticTetra()
                elif cellType == 'C3D15':
                    cell = vtk.vtkQuadraticWedge()
                elif cellType == 'C3D20' or cellType == 'C3D20R':
                    cell = vtk.vtkQuadraticHexahedron()
                elif cellType == 'C2D3':
                    cell = vtk.vtkTriangle()
                elif cellType == 'C2D4':
                    cell = vtk.vtkQuad()
                elif cellType == 'C2D6':
                    cell = vtk.vtkQuadraticTriangle()
                elif cellType == 'C2D8':
                    cell = vtk.vtkQuadraticQuad()
                else:
                    raise ValueError(f'Unknown cell type: {cellType}')
                for j, nodeId in enumerate(nodes):
                    cell.GetPointIds().SetId(j, nodeId)

                self.ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
            else:
                self.ugrid.InsertNextCell(0, vtk.vtkIdList())
    def _write_vtk(self):
        # Write the geometry and solution data to a VTK file
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(self.output_filename)
        writer.SetInputData(self.ugrid)

        # Add solution fields to the VTK file
        fields_dict = self.solution.list_fields_reduced()

        for field_name, field_func in fields_dict.items():
            field_data = field_func()

            # Check if data is node-based or element-based
            if len(field_data) == len(self.geometry.nodes):  # Node-based data
                vtk_data = numpy_to_vtk(field_data)
                vtk_data.SetName(field_name)
                self.ugrid.GetPointData().AddArray(vtk_data)
            elif len(field_data) == len(self.geometry.elements):  # Element-based data
                vtk_data = numpy_to_vtk(field_data)
                vtk_data.SetName(field_name)
                self.ugrid.GetCellData().AddArray(vtk_data)
            else:
                print(f"Warning: Field '{field_name}' size does not match nodes or elements.")

        writer.Write()

    def convert(self):
        # Generate unstructured grid from geometry
        self._generate_ugrid()
        # Write the unstructured grid with solution fields to VTK file
        self._write_vtk()


def main():
    parser = argparse.ArgumentParser(description="Convert geometry and solution data to VTK format.")
    parser.add_argument("geometry_path", help="Path to the geometry input file (.inp).")
    parser.add_argument(
        "--solution_path",
        help="Path to the solution file (.res). If not provided, deduced from the input file by replacing '.inp' with '.res'.",
        default=None,
    )
    parser.add_argument(
        "--output",
        help="Output VTK file name. If not provided, deduced from the input file by replacing '.inp' with '.vtk'.",
        default=None
    )

    args = parser.parse_args()

    # Create the converter instance and perform the conversion
    converter = Convert(geometry_path=args.geometry_path, solution_path=args.solution_path, output_filename=args.output)
    converter.convert()


if __name__ == '__main__':
    main()
