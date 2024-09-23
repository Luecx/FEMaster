import vtk
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk
from .geometry import Geometry
from .solution import Solution

class Convert:
    def __init__(self, geometry_path, solution_path, output_filename="output.vtk"):
        self.geometry = Geometry.read_input_deck(geometry_path)
        self.solution = Solution.open(solution_path)
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

if __name__ == '__main__':
    # Example usage:
    geometry_path = "res/haken_analysis_2/v1.inp"
    solution_path = "res/haken_analysis_2/v1.inp.res"
    output_filename = "output1.vtk"

    converter = Convert(geometry_path, solution_path, output_filename)
    converter.convert()
