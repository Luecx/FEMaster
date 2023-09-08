import geometry
import vtk
import numpy as np
import matplotlib.pyplot as plt
import solution
from vtk.util.numpy_support import numpy_to_vtk

class Viewer:
    def __init__(self):
        self.geometry = None
        self.data = None
        self.data_type = None
        self.data_min = None
        self.data_max = None
        self.color_scheme = None

        self.displacement = None
        self.animate = False

        self.cos = False

        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.renderWindowInteractor.SetRenderWindow(self.renderWindow)
        self.actor = vtk.vtkActor()
        self.renderer.AddActor(self.actor)

        self.renderer.SetBackground(0.9,0.9,0.9)

        # self.renderer.GradientBackgroundOn()  # Enable gradient background
        # self.renderer.SetBackground(0.6, 0.6, 0.66)  # Bottom color
        # self.renderer.SetBackground2(0.5, 0.5, 0.5)  # Top color


    def set_geometry(self, geometry):
        self.geometry = geometry

    def set_data(self, type='node', data=None):
        if data is not None:
            self.data = numpy_to_vtk(data)
        else:
            self.data = None
        self.data_type = type
        if data is not None:
            self.set_data_range(np.min(data), np.max(data))

    def set_data_range(self, min, max, percentile=False):
        if percentile and self.data is not None:
            self.data_min = np.percentile(self.data, min * 100)
            self.data_max = np.percentile(self.data, max * 100)
        else:
            self.data_min = min
            self.data_max = max

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

    def coordinate_system(self):
        self.cos = True

    def animate(self):
        self.animate = True

    def _update_geometry(self, obj, event):
        # Compute the displacement at the current time
        current_time = obj.GetTimerDuration()
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


    def _visualize_geom(self):
        # Create the unstructured grid object
        self.ugrid = vtk.vtkUnstructuredGrid()

        # Set the points
        points = vtk.vtkPoints()
        for i, node in enumerate(self.geometry.nodes):
            if node:
                if not self.animate and self.displacement is not None:
                    points.InsertNextPoint(node + self.displacement[i])
                else:
                    points.InsertNextPoint(node)
            else:
                points.InsertNextPoint([0,0,0])
        self.ugrid.SetPoints(points)

        # Set the cells
        for element in self.geometry.elements:
            if element:
                cellType = element['type']
                nodes = element['nodes']
                if cellType in ['C3D4', 'C3D10']:
                    cell = vtk.vtkTetra()
                    num_nodes = 4
                elif cellType in ['C3D6', 'C3D15']:
                    cell = vtk.vtkWedge()
                    num_nodes = 6
                elif cellType in ['C3D8', 'C3D20']:
                    cell = vtk.vtkHexahedron()
                    num_nodes = 8
                else:
                    raise ValueError(f'Unknown cell type: {cellType}')
                for i, nodeId in enumerate(nodes[:num_nodes]):
                    cell.GetPointIds().SetId(i, nodeId)
                self.ugrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

        # Map the unstructured grid to colors
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(self.ugrid)
        mapper.SetInterpolateScalarsBeforeMapping(1)
        self.actor.SetMapper(mapper)

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

    def _visualize_cos(self):
        self.axes = vtk.vtkAxesActor()
        self.renderer.AddActor(self.axes)

    def mainloop(self):
        if self.geometry:
            self._visualize_geom()

        if self.cos:
            self._visualize_cos()

        # Start the render loop
        self.renderWindow.Render()

        # Create a timer event
        if self.displacement is not None and self.animate:
            self.renderWindowInteractor.AddObserver('TimerEvent', self._update_geometry)
            self.renderWindowInteractor.CreateRepeatingTimer(100)

        self.renderWindowInteractor.Start()



def mises(stress):
    sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
    mises = np.sqrt(0.5 * ((sigma_x - sigma_y) ** 2 + (sigma_y - sigma_z) ** 2 +
                           (sigma_z - sigma_x) ** 2 + 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_zx ** 2)))
    return mises


def principal(stress):
    sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
    num_points = sigma_x.shape[0]
    principal_stresses = np.zeros((num_points, 3))

    for i in range(num_points):
        stress_tensor = np.array([[sigma_x[i], tau_xy[i], tau_zx[i]],
                                  [tau_xy[i], sigma_y[i], tau_yz[i]],
                                  [tau_zx[i], tau_yz[i], sigma_z[i]]])
        eigvals, _ = np.linalg.eig(stress_tensor)
        principal_stresses[i] = np.sort(eigvals)

    return principal_stresses


def signed_mises(stress):
    mises_stress = mises(stress)
    principal_stresses_val = principal(stress)
    max_absolute_principal_stress = np.max(np.abs(principal_stresses_val), axis=1)
    sign = np.sign(np.sum(principal_stresses_val * (np.abs(principal_stresses_val) == max_absolute_principal_stress[:, None]), axis=1))
    signed_mises_stress = sign * mises_stress
    return signed_mises_stress



# Create a model

geom = geometry.Geometry.read_input_deck("../../topo/runs/topo_20230907_08-45-11/redesign/opt.inp")
solution = solution.Solution.open("../../topo/runs/topo_20230907_08-45-11/redesign/opt.inp.res")

disp = solution.loadcases["1"]["DISPLACEMENT"]
strain = solution.loadcases["1"]["STRAIN"]
stress = solution.loadcases["1"]["STRESS"]
mises = signed_mises(stress)


# viewer
viewer = Viewer()
viewer.set_geometry(geom)
viewer.set_data(type='node', data=mises)
viewer.set_data_range(0.01, 0.99, percentile=True)
viewer.set_colorscheme('jet')
viewer.set_displacement(disp[:,:3] * 0.1)
viewer.coordinate_system()
viewer.mainloop()