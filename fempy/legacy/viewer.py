
# import geometry and solution in the same folder
from .geometry import Geometry
from .solution import Solution

import vtk
import numpy as np
import math
import matplotlib.pyplot as plt
import time
import sys
from vtk.util.numpy_support import numpy_to_vtk
from vtk import *
import argparse

class CustomInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, parent=None):
        self.AddObserver("MiddleButtonPressEvent", self.middle_button_press_event)
        self.AddObserver("MiddleButtonReleaseEvent", self.middle_button_release_event)

    def middle_button_press_event(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # Get the world coordinates of the picked point
        pickPosition = picker.GetPickPosition()

        if picker.GetActor():
            # Set the new focal point
            self.GetDefaultRenderer().GetActiveCamera().SetFocalPoint(pickPosition)
            self.GetDefaultRenderer().ResetCameraClippingRange()

        self.OnMiddleButtonDown()
        return

    def middle_button_release_event(self, obj, event):
        self.OnMiddleButtonUp()
        return



class Viewer:
    def __init__(self):

        # data visual of main body
        self.geometry = None
        self.elem_mask = None
        self.data = None
        self.data_type = None
        self.data_min = None
        self.data_max = None
        self.color_scheme = None
        self.model_size = None
        self.boundaries = False

        # animation and display of displacements
        self.displacement = None
        self.animate = False
        self.start_time = time.time()

        # true if drawing a cos
        self.cos = 0

        # grid spacing
        self.grid_xy = 0

        # vtk related stuff
        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.renderWindowInteractor.SetRenderWindow(self.renderWindow)

        self.customInteractorStyle = CustomInteractorStyle()
        self.customInteractorStyle.SetDefaultRenderer(self.renderer)
        self.renderWindowInteractor.SetInteractorStyle(self.customInteractorStyle)

        self.actor = vtk.vtkActor()
        self.renderer.AddActor(self.actor)

        self.renderer.SetBackground(0.9,0.9,0.9)

        # self.renderer.GradientBackgroundOn()  # Enable gradient background
        # self.renderer.SetBackground(0.6, 0.6, 0.66)  # Bottom color
        # self.renderer.SetBackground2(0.5, 0.5, 0.5)  # Top color


    def set_geometry(self, geometry):
        self.geometry = geometry
        filtered_nodes = [node for node in geometry.nodes if node is not None]

        print("setting geometry with ", len(filtered_nodes), " nodes")

        # Extract the x, y, and z coordinates
        x_coords = [node[0] for node in filtered_nodes]
        y_coords = [node[1] for node in filtered_nodes]
        z_coords = [node[2] for node in filtered_nodes]

        # Compute the extents
        x_extent = max(x_coords) - min(x_coords)
        y_extent = max(y_coords) - min(y_coords)
        z_extent = max(z_coords) - min(z_coords)

        # Find the largest extent
        self.model_size = max(x_extent, y_extent, z_extent)

    def set_element_mask(self, mask):
        self.elem_mask = mask

    def set_data(self, type='node', data=None):
        if data is not None:
            self.data = numpy_to_vtk(data)
        else:
            self.data = None
        self.data_type = type
        if data is not None:
            self.set_data_range(np.min(data), np.max(data))

    def set_boundaries(self, value=True):
        self.boundaries = value

    def set_data_range(self, min, max, percentile=False):
        if percentile and self.data is not None:
            self.data_min = np.percentile(self.data, min * 100)
            self.data_max = np.percentile(self.data, max * 100)
        else:
            self.data_min = min
            self.data_max = max

    def set_grid_xy(self, spacing=None):
        if spacing is not None:
            self.grid_xy = spacing
        elif self.model_size is not None:
            self.grid_xy = 10 ** math.floor(math.log10(self.model_size / 3))

    def set_colorscheme(self, scheme):
        self.color_scheme = scheme
        lut = vtk.vtkLookupTable()

        # Create a color map with 20 colors
        num_colors = 20
        lut.SetNumberOfTableValues(num_colors)
        cmap = plt.get_cmap(scheme)
        for i in range(num_colors):
            # Set the RGB values for each entry in the table
            r, g, b, _ = cmap(i / (num_colors - 1))
            lut.SetTableValue(i, r, g, b, 1.0)

        lut.SetNumberOfTableValues(num_colors)
        lut.Build()
        self.lut = lut

    def set_displacement(self, displacement):
        self.displacement = displacement

    def coordinate_system(self, size=None):
        if size is not None:
            self.cos = size
        elif self.model_size is not None:
            self.cos = 10 ** math.floor(math.log10(self.model_size / 3))

    def add_animation(self):
        self.animate = True

    def _update_geometry(self, obj, event):
        # Compute the displacement at the current time
        current_time = time.time() - self.start_time
        current_disp = self.displacement * np.sin(current_time)
        # Update the geometry
        points = vtk.vtkPoints()
        for node, disp in zip(self.geometry.nodes, current_disp):
            if node:
                points.InsertNextPoint(node + disp)
            else:
                points.InsertNextPoint([0, 0, 0])
        self.ugrid.SetPoints(points)
        # Render the scene
        self.renderWindow.Render()

    def write_to_stl(self, filename="output.stl"):
        if not hasattr(self, 'ugrid'):
            self._generate_ugrid()

        # Convert vtkUnstructuredGrid to vtkPolyData
        geometryFilter = vtk.vtkGeometryFilter()
        geometryFilter.SetInputData(self.ugrid)
        geometryFilter.Update()

        polyData = geometryFilter.GetOutput()

        # Create an STL writer and set the filename
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(filename)
        stlWriter.SetInputData(polyData)
        stlWriter.Write()


    def _generate_ugrid(self):
        # Create the unstructured grid object
        self.ugrid = vtk.vtkUnstructuredGrid()
        # Set the points
        points = vtk.vtkPoints()
        for i, node in enumerate(self.geometry.nodes):
            if node:
                if not self.animate and self.displacement is not None:
                    points.InsertNextPoint(node[0] + self.displacement[i][0], node[1] + self.displacement[i][1], node[2] + self.displacement[i][2])
                else:
                    points.InsertNextPoint(node)
            else:
                points.InsertNextPoint([0, 0, 0])
        self.ugrid.SetPoints(points)

        # Set the cells
        for i, element in enumerate(self.geometry.elements):
            if element and (self.elem_mask is None or self.elem_mask[i]):
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
                cell = vtk.vtkTetra()
                cell.GetPointIds().SetId(0, 0)
                cell.GetPointIds().SetId(1, 0)
                cell.GetPointIds().SetId(2, 0)
                cell.GetPointIds().SetId(3, 0)
                self.ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

    def _visualize_geom(self):
        # Create the unstructured grid object
        if not hasattr(self, 'ugrid'):
            self._generate_ugrid()

        # Map the unstructured grid to colors
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(self.ugrid)
        mapper.SetInterpolateScalarsBeforeMapping(1)

        # Set up the main actor with the mapper
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetRepresentationToSurface()  # Set the main representation to surface

        # If data, data_min, data_max, and color_scheme are provided, map it to colors and display color bar
        if self.data is not None and self.data_min is not None and self.data_max is not None and self.color_scheme is not None:
            # Create a lookup table to map point data to colors
            self.lut.SetRange(self.data_min, self.data_max)
            mapper.SetLookupTable(self.lut)

            if self.data_type == 'node':
                self.ugrid.GetPointData().SetScalars(self.data)
                mapper.SetScalarRange(self.data_min, self.data_max)
            elif self.data_type == 'element':
                self.ugrid.GetCellData().SetScalars(self.data)
                mapper.SetScalarRange(self.data_min, self.data_max)
            mapper.SetLookupTable(self.lut)

            # Show color bar
            colorbar = vtk.vtkScalarBarActor()
            colorbar.SetLookupTable(self.lut)
            colorbar.SetTitle("Color Map")
            colorbar.SetNumberOfLabels(20)
            self.renderer.AddActor(colorbar)

        # Add the wireframe overlay if boundaries are enabled
        if self.boundaries:
            wireframe_mapper = vtk.vtkDataSetMapper()
            wireframe_mapper.SetInputData(self.ugrid)
            wireframe_actor = vtk.vtkActor()
            wireframe_actor.SetMapper(wireframe_mapper)
            wireframe_actor.GetProperty().SetRepresentationToWireframe()
            wireframe_actor.GetProperty().SetColor(0, 0, 0)  # Set color to black
            wireframe_actor.GetProperty().SetLineWidth(2.0)  # Set line width

            # Add the wireframe actor to the renderer
            self.renderer.AddActor(wireframe_actor)

        # Ensure background and other settings are set if not already
        self.renderer.SetBackground(0.9, 0.9, 0.9)  # This can be changed to any color

    def _visualize_cos(self):
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(self.cos, self.cos, self.cos)
        self.renderer.AddActor(self.axes)

    def _visualize_grid(self):
        # create grid lines
        for i in range(-5, 6):
            line_x = vtk.vtkLineSource()
            line_x.SetPoint1(i*self.grid_xy, -5*self.grid_xy, 0)
            line_x.SetPoint2(i*self.grid_xy, 5*self.grid_xy, 0)

            line_y = vtk.vtkLineSource()
            line_y.SetPoint1(-5*self.grid_xy, i*self.grid_xy, 0)
            line_y.SetPoint2(5*self.grid_xy, i*self.grid_xy, 0)

            mapper_x = vtk.vtkPolyDataMapper()
            mapper_x.SetInputConnection(line_x.GetOutputPort())
            actor_x = vtk.vtkActor()
            actor_x.SetMapper(mapper_x)

            mapper_y = vtk.vtkPolyDataMapper()
            mapper_y.SetInputConnection(line_y.GetOutputPort())
            actor_y = vtk.vtkActor()
            actor_y.SetMapper(mapper_y)

            self.renderer.AddActor(actor_x)
            self.renderer.AddActor(actor_y)


    def mainloop(self):
        if self.geometry:
            self._visualize_geom()

        if self.grid_xy:
            self._visualize_grid()

        if self.cos:
            self._visualize_cos()

        # Start the render loop
        self.renderWindow.Render()

        # Create a timer event
        if self.displacement is not None and self.animate:
            self.renderWindowInteractor.AddObserver('TimerEvent', self._update_geometry)
            self.renderWindowInteractor.CreateRepeatingTimer(40)

        self.renderWindowInteractor.Start()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Visualize geometry and solutions")
    parser.add_argument("--geometry", type=str, help="Path to the geometry input deck")

    parser.add_argument("--solution", type=str, help="Path to the solution file")
    parser.add_argument("--loadcase", type=str, default='1', help="Specify loadcase")
    parser.add_argument("--field", type=str, help="Field to display")

    parser.add_argument("--disp_scalar", type=float, default=1.0, help="Displacement scaling factor")
    parser.add_argument("--disp_field", type=str, default="displacement_xyz", help="Field to use for displacement")
    parser.add_argument("--displacement", action='store_true', help="If set, enable displacements")

    parser.add_argument("--boundaries", action='store_true', help="If set, enabled boundaries of elements being displayed")

    parser.add_argument("--datarange", nargs=3, help="Data range (min, max, percentile). Use True/False for percentile.")

    args = parser.parse_args()

    viewer = Viewer()

    if args.geometry:
        viewer.set_geometry(Geometry.read_input_deck(args.geometry))
        viewer.set_boundaries(args.boundaries)

    if args.solution:
        sol = Solution.open(args.solution)
        available_fields = sol.list_fields(loadcase=args.loadcase)

        if args.field:
            if args.field not in available_fields:
                print(f"Error: The field '{args.field}' is not available for the given loadcase.")
                print("Possible fields:")
                for k in available_fields:
                    print(f"\t{k}")
                print("Raw fields:")
                print(sol)
                sys.exit(1)

            print(f"Displaying field: {args.field}")
            print("geometry nodes: ", len(viewer.geometry.nodes))
            print("values in field: ", len(available_fields[args.field]()))

            if len(viewer.geometry.nodes) == len(available_fields[args.field]()):
                viewer.set_data(type='node', data=available_fields[args.field]())
            else:
                viewer.set_data(type='element', data=available_fields[args.field]())
            viewer.set_colorscheme('jet')

        if args.datarange:
            min_val, max_val, percentile = args.datarange
            viewer.set_data_range(float(min_val), float(max_val), percentile=(percentile.lower() == 'true'))

        if args.displacement:
            # Use the specified displacement field or default to 'displacement_xyz'
            disp_field = args.disp_field
            if disp_field not in available_fields:
                print(f"Error: The field '{disp_field}' is not available for the given loadcase.")
                print("Possible fields:")
                for k in available_fields:
                    print(f"\t{k}")
                sys.exit(1)
            # Check that the displacement field has exactly 3 columns
            disp_data = available_fields[disp_field]()
            if disp_data.shape[1] != 3:
                print(f"Error: The displacement field '{disp_field}' must have exactly 3 columns.")
                sys.exit(1)
            # Apply the displacement with the specified scalar factor
            viewer.set_displacement(disp_data * args.disp_scalar)



    viewer.coordinate_system()
    viewer.set_grid_xy()
    viewer.mainloop()
