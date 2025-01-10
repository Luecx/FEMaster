
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
        self.curvature = None

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
    def __init__(self, points, type="outer", clockwise=False):
        self.points = points
        self.type = type
        self.clockwise = clockwise
        self.compute_geometry_data()

    def compute_geometry_data(self):
        for i in range(len(self.points)):
            p_prev = self.points[(i - 1 + len(self.points)) % len(self.points)]
            p = self.points[i]
            p_next = self.points[(i + 1) % len(self.points)]

            # Convert to numpy arrays
            _p_prev = np.array([p_prev.x, p_prev.y])
            _p = np.array([p.x, p.y])
            _p_next = np.array([p_next.x, p_next.y])

            # Compute the normal
            tangent = _p_next - _p_prev
            normal = np.array([tangent[1], -tangent[0]])

            # if the boundary is clockwise, invert the normal
            if self.clockwise:
                normal = -normal
            # if the boundary is inner, invert the normal
            if self.type == "inner":
                normal = -normal

            normal = normal / np.linalg.norm(normal)

            # Compute the curvature
            area = 0.5 * np.abs(np.cross(_p_prev - _p, _p_next - _p))
            side1 = np.linalg.norm(_p_prev - _p)
            side2 = np.linalg.norm(_p_next - _p)
            side3 = np.linalg.norm(_p_next - _p_prev)
            curvature = 4 * area / (side1 * side2 * side3)

            p.normal = normal
            p.curvature = curvature

    def update_geometry(self, movements):
        for point in self.points:
            if point.rigid:
                continue
            point.x += movements[point.id] * point.normal[0]
            point.y += movements[point.id] * point.normal[1]
        self.compute_geometry_data()

    def redistribute_nodes(self):
        # Keep a copy of the old positions
        old_positions = [(p.x, p.y) for p in self.points]

        # Compute the total length of the structure
        total_length = 0
        distances = []
        for i in range(len(self.points)):
            p1 = np.array([self.points[i].x, self.points[i].y])
            p2 = np.array([self.points[(i + 1) % len(self.points)].x, self.points[(i + 1) % len(self.points)].y])
            distance = np.linalg.norm(p2 - p1)
            distances.append(distance)
            total_length += distance

        # Compute the target spacing
        num_points = len(self.points)
        target_spacing = total_length / num_points

        # Redistribute the nodes
        new_positions = [old_positions[0]]  # Start with the first node
        accumulated_distance = 0
        current_index = 0

        for i in range(1, num_points):
            # Target distance to the next node
            target_distance = i * target_spacing

            while accumulated_distance + distances[current_index] < target_distance:
                accumulated_distance += distances[current_index]
                current_index = (current_index + 1) % len(self.points)

            # Interpolate between the current and next node
            p1 = np.array(old_positions[current_index])
            p2 = np.array(old_positions[(current_index + 1) % len(old_positions)])
            segment_length = distances[current_index]
            remaining_distance = target_distance - accumulated_distance
            t = remaining_distance / segment_length
            new_position = p1 + t * (p2 - p1)
            new_positions.append((new_position[0], new_position[1]))

        # Update the points with the new positions
        for i, point in enumerate(self.points):
            point.x, point.y = new_positions[i]

        # Recompute geometry data
        self.compute_geometry_data()

        return old_positions, new_positions

class Problem2D:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self._numerate_points()
        self._compute_geometric_data()

    def _create_geom_2d(self):
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
        return geom_2d

    def _create_geom_3d(self):

        geom_2d = self._create_geom_2d()
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

    def _numerate_points(self):
        global_id = 0
        for boundary in self.boundaries:
            for _, point in enumerate(boundary.points):
                point.id = global_id
                global_id += 1

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

    def _compute_geometric_data(self):
        for boundary in self.boundaries:
            boundary.compute_geometry_data()

    def _run(self, geom_file):
        # write a .inp file, get time in ns
        timehash = time.time_ns()

        try:
            geom_file.write_input_deck(f"{timehash}.inp")
            self._attach_bc(f"{timehash}.inp")

            # run the solver and wait
            sp.run(["../../bin/FEMaster", f"{timehash}.inp"], stdout=sp.DEVNULL, stderr=sp.DEVNULL)

            # read the results
            results = Solution.open(f"{timehash}.res", loadingbar=False)
            # read the mises stress
            mises = results.list_fields()['mises']()
            self._store_mises(mises)

            os.remove(f"{timehash}.inp")
            os.remove(f"{timehash}.res")
        except Exception as e:
            pass
            os.remove(f"{timehash}.inp")
            os.remove(f"{timehash}.res")

    def _store_mises(self, mises):
        for boundary in self.boundaries:
            for point in boundary.points:
                point.fem_mises = 0
                id = point.fem_ids[0]
                point.fem_mises = mises[id]

    def run(self):
        geom_3d = self._create_geom_3d()
        self._run(geom_3d)

    def update_geom(self, movements):
        for boundary in self.boundaries:
            boundary.update_geometry(movements)
            boundary.redistribute_nodes()

    def max_mises(self):
        try:
            return max([point.fem_mises for point in self])
        except:
            return 0
    def get_inputs(self):
        inputs = []
        max_mises = self.max_mises()
        for boundary in self.boundaries:
            for point in boundary.points:
                inputs.append([point.normal[0], point.normal[1], point.curvature, point.fem_mises / max_mises])
        return np.array(inputs)

    # iterator over all points
    def __iter__(self):
        for boundary in self.boundaries:
            for point in boundary.points:
                yield point

    def redistribute_nodes(self):
        for boundary in self.boundaries:
            boundary.redistribute_nodes()

    def plot(self):
        # plot the nodes with their normals
        import matplotlib.pyplot as plt
        for boundary in self.boundaries:
            for point in boundary.points:
                plt.plot(point.x, point.y, 'ro')
                plt.quiver(point.x, point.y, point.normal[0], point.normal[1])
        plt.show()

    def adjacency_matrix(self):
        # Initialize an adjacency matrix
        num_points = len([point for point in self])
        adjacency = np.zeros((num_points, num_points))

        # Fill the adjacency matrix
        for boundary in self.boundaries:
            for localid, point in enumerate(boundary.points):
                localid_prev = (localid - 1 + len(boundary.points)) % len(boundary.points)
                localid_next = (localid + 1) % len(boundary.points)

                id = point.id
                id_prev = boundary.points[localid_prev].id
                id_next = boundary.points[localid_next].id

                adjacency[id, id_prev] = 1
                adjacency[id, id_next] = 1
                adjacency[id_prev, id] = 1
                adjacency[id_next, id] = 1
                adjacency[id, id] = 1

        return adjacency


def create_geom():

    # Create the geometry
    outer_points = []
    inner_points = []
    inner_points2 = []

    for i in range(10):
        outer_points.append(Point2D(id=i, x = i, y = 0, rigid=True))
    for i in range(10):
        outer_points.append(Point2D(id=i+10, x = 10, y = i, rigid=True))
    for i in range(10):
        outer_points.append(Point2D(id=i+20, x = 10-i, y = 10, rigid=True))
    for i in range(10):
        outer_points.append(Point2D(id=i+30, x = 0, y = 10-i, rigid=True))

    for i in range(11):
        # support in x and y
        outer_points[i].bc(disp_y = 0)
        outer_points[i+10].bc(disp_x = 0)
        outer_points[i+20].bc(force_x=0, force_y=1)
        outer_points[(i+30)%40].bc(force_x=-1, force_y=0)

    N_inner = 10
    # inner points is simply a circle
    for i in range(N_inner):
        inner_points.append(Point2D(id=i+40, x = 5 + 3*np.cos(i/N_inner*2*np.pi), y = 5 + 3*np.sin(i/N_inner*2*np.pi)))
    # for i in range(N_inner):
    #     inner_points2.append(Point2D(id=i+60, x = 5 + 1.5*np.cos(i/N_inner*2*np.pi), y = 3 + 1.5*np.sin(i/N_inner*2*np.pi)))

    outer_boundary = Boundary2D(outer_points, type="outer")
    inner_boundary = Boundary2D(inner_points, type="inner")
    # inner_boundary2 = Boundary2D(inner_points2, type="inner")
    problem = Problem2D([outer_boundary, inner_boundary])


    return problem

# problem.run()
# print(problem.max_mises())
# problem.plot()
#
# for i in range(5):
#     inputs = problem.get_inputs()
#     # try to move as a function of the mises, if mises > 0.5, move in the normal direction
#     problem.update_geom({point.id: 0.5*(inputs[point.id][3] - 0.5) for point in problem})
#     problem.run()
#     geom2d = problem._create_geom_2d()
#     geom2d.plot_2d()
#     print(problem.max_mises())
