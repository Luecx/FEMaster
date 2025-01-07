
from fempy.geometry import *
from fempy.solution import *

import numpy as np
import subprocess as sp

import time
import os


class Point2D:
    def __init__(self, id, x, y, rigid=False):
        self.id = id
        self.x = x
        self.y = y
        self.rigid = rigid

        # boundary condition
        self.bc_disp_x = None
        self.bc_disp_y = None
        self.bc_force_x = None
        self.bc_force_y = None

        # results
        self.fem_ids = None
        self.fem_mises = None
        self.normal = None

    def bc(self, disp_x=None, disp_y=None, force_x=None, force_y=None):
        self.bc_disp_x = disp_x
        self.bc_disp_y = disp_y
        self.bc_force_x = force_x
        self.bc_force_y = force_y

    def __str__(self):
        s = "Point2D(" + str(self.id) + ", " + str(self.x) + ", " + str(self.y) + ")"
        s += "  FEM ID: " + str(self.fem_ids)
        s += "  FEM Mises: " + str(self.fem_mises)
        return s

class Boundary2D:
    def __init__(self, points):
        self.points = points

class Problem2D:
    def __init__(self, boundaries):
        self.boundaries = boundaries

    def _create_geom(self):
        segment_groups = []
        # for each boundary, create a segment group with segments
        for boundary in self.boundaries:
            segments = []
            for i in range(len(boundary.points)):
                p1 = boundary.points[i]
                p2 = boundary.points[(i+1) % len(boundary.points)]
                p1xy = np.array([p1.x, p1.y])
                p2xy = np.array([p2.x, p2.y])
                segments.append(StraightSegment(p1xy, p2xy, n_subdivisions=1))
            segment_group = SegmentGroup(segments)
            segment_groups.append(segment_group)

        geom_2d = Geometry.mesh_interior(segment_groups)
        geom_3d = geom_2d.extruded(1, 1)

        # match the ids from the geometry to each node
        for boundary in self.boundaries:
            for point in boundary.points:
                point.fem_ids = []

                for id, node in enumerate(geom_3d.nodes):
                    if node is None:
                        continue

                    if np.allclose([point.x, point.y, 0], node):
                        point.fem_ids.append(id)

        return geom_3d.as_second_order()

    def _attach_bc(self, geom_file):
        with open(geom_file, "a") as file:

            # write support
            file.write("*SUPPORT, SUPPORT_COLLECTOR=SUPPS\n")
            for boundary in self.boundaries:
                for point in boundary.points:
                    for id in point.fem_ids:
                        if point.bc_disp_x is not None or point.bc_disp_y is not None:
                            file.write(f"{id}, "
                                       f"{point.bc_disp_x if point.bc_disp_x is not None else ''}, "
                                       f"{point.bc_disp_y if point.bc_disp_y is not None else ''}\n")
            file.write("*CLOAD, LOAD_COLLECTOR=LOADS\n")
            for boundary in self.boundaries:
                for point in boundary.points:
                    for id in point.fem_ids:
                        if point.bc_force_x is not None or point.bc_force_y is not None:
                            file.write(f"{id}, "
                                       f"{point.bc_force_x if point.bc_force_x is not None else ''}, "
                                       f"{point.bc_force_y if point.bc_force_y is not None else ''}\n")

            # write material
            file.write("*MATERIAL, NAME=MAT1\n")
            file.write("*ELASTIC, TYPE=ISOTROPIC\n")
            file.write("200000., 0.3\n")
            # solid section
            file.write("*SOLID SECTION, ELSET=EALL, MATERIAL=MAT1\n")
            # loadcase
            file.write("*LOAD CASE, type=linear static\n")
            file.write("*LOAD\n")
            file.write("LOADS\n")
            file.write("*SUPPORT\n")
            file.write("SUPPS\n")
            file.write("*END")

    def _run(self, geom_file):
        # write a .inp file, get time in ns
        timehash = time.time_ns()

        try:
            geom_file.write_input_deck(f"{timehash}.inp")
            self._attach_bc(f"{timehash}.inp")

            # run the solver and wait
            sp.run(["../../bin/FEMaster", f"{timehash}.inp"])

            # read the results
            results = Solution.open(f"{timehash}.res", loadingbar=False)
            # read the mises stress
            mises = results.list_fields()['mises']()
            self._store_mises(mises)

            # os.remove(f"{timehash}.inp")
            # os.remove(f"{timehash}.res")
        except Exception as e:
            pass
            # os.remove(f"{timehash}.inp")
            # os.remove(f"{timehash}.res")

    def _store_mises(self, mises):
        for boundary in self.boundaries:
            for point in boundary.points:
                point.fem_mises = 0
                id = point.fem_ids[0]
                point.fem_mises = mises[id]

# Create the geometry
outer_points = []
inner_points = []

for i in range(10):
    outer_points.append(Point2D(id=i, x = i, y = 0))
for i in range(10):
    outer_points.append(Point2D(id=i+10, x = 10, y = i))
for i in range(10):
    outer_points.append(Point2D(id=i+20, x = 10-i, y = 10))
for i in range(10):
    outer_points.append(Point2D(id=i+30, x = 0, y = 10-i))

for i in range(11):
    # support in x and y
    outer_points[i].bc(disp_x = 0, disp_y = 0)
    outer_points[i+20].bc(force_x=0, force_y=-1)

# inner points is simply a circle
for i in range(20):
    inner_points.append(Point2D(id=i+40, x = 5 + 2*np.cos(i/20*2*np.pi), y = 5 + 2*np.sin(i/20*2*np.pi)))


outer_boundary = Boundary2D(outer_points)
inner_boundary = Boundary2D(inner_points)
problem = Problem2D([outer_boundary, inner_boundary])

geom_3d = problem._create_geom()
problem._run(geom_3d)

# plot the points, each with their mises stress as a color
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
for boundary in problem.boundaries:
    for point in boundary.points:
        plt.plot(point.x, point.y, 'o', color=plt.cm.viridis(point.fem_mises/10))
plt.show()