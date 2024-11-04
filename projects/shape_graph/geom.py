import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline

# add "../../python/" to the path
import sys
import os
sys.path.append("../../")

from fempy.geometry import *
from fempy.solution import *



class Point:
    def __init__(self, pos, tangent=None, curvature=None, stress=None, id=None):
        self.pos = pos
        self.tangent = tangent
        self.curvature = curvature
        self.stress = stress
        self.id = None

    def compute_curvature(self, before, after):
        # Calculate the curvature of a point
        # given the previous, current, and next points
        v1 = self.pos - before.pos
        v2 = self.pos - after.pos
        v3 = after.pos - before.pos

        # extend to 3d
        v1 = np.append(v1, 0)
        v2 = np.append(v2, 0)
        v3 = np.append(v3, 0)

        area = np.linalg.norm(np.cross(v1, v2)) / 2
        curv = 4 * area / np.linalg.norm(v1) / np.linalg.norm(v2) / np.linalg.norm(v3)
        return curv

    def rotated(self, angle):
        return Point(
            pos=np.array([
                self.pos[0] * np.cos(angle) - self.pos[1] * np.sin(angle),
                self.pos[0] * np.sin(angle) + self.pos[1] * np.cos(angle)
            ]),
            tangent=np.array([
                self.tangent[0] * np.cos(angle) - self.tangent[1] * np.sin(angle),
                self.tangent[0] * np.sin(angle) + self.tangent[1] * np.cos(angle)
            ]),
            curvature=self.curvature,
            stress=self.stress,
            id=self.id
        )

    def __repr__(self):
        return f"Point" \
               f"\n    position  = {self.pos}, " \
               f"\n    tangent   = {self.tangent}, " \
               f"\n    curvature = {self.curvature}, " \
               f"\n    stress    = {self.stress}, " \
               f"\n    id        = {self.id})"

class PointGroup:
    def __init__(self, points):
        self.points = points

    def compute_curvature(self):
        for i in range(1, len(self.points)-1):
            p0 = self.points[i-1]
            p1 = self.points[i]
            p2 = self.points[i+1]
            self.points[i].curvature = p1.compute_curvature(p0, p2)

    def compute_tangents(self):
        for i in range(1, len(self.points)-1):
            p0 = self.points[i-1]
            p1 = self.points[i]
            p2 = self.points[i+1]
            v1 = p1.pos - p0.pos
            v2 = p2.pos - p1.pos
            self.points[i].tangent = (v1 + v2) / np.linalg.norm(v1 + v2)

    def compute_normals(self):
        # towards the center of curvature
        for i in range(1, len(self.points)-1):
            p0 = self.points[i-1]
            p1 = self.points[i]
            p2 = self.points[i+1]

            v0 = p0.pos - p1.pos
            v2 = p2.pos - p1.pos

            normal = v0 + v2
            normal /= np.linalg.norm(normal)



class Model:

    def __init__(self, point1=None, point2=None, point3=None, num_points=35,id=0):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3

        self.id = id

        self.num_points  = num_points
        self.base_points = None
        self.point_group = None

        self._init_default_points()
        self._create_point_group()



    def _init_default_points(self):
        np.random.seed(self.id)

        def _rand(center, deviation):
            return center + np.random.uniform(-deviation, deviation)

        def _pos(x_d, y_d):
            return np.array([_rand(*x_d), _rand(*y_d)])

        def _rot(vec, angle):
            return np.array([
                vec[0] * np.cos(angle) - vec[1] * np.sin(angle),
                vec[0] * np.sin(angle) + vec[1] * np.cos(angle)
            ])

        #     p2_x_target, p2_y_target = 1, 0.5
        #     p2_x_deviation, p2_y_deviation = 1, 0.1
        #
        #     p3_x_target, p3_y_target = 0, 1
        #     p3_x_deviation, p3_y_deviation = 0.5, 0.1
        #
        #     # Angle deviation controls
        #     t1_angle_dev = 0.1  # Angle deviation for P1 tangent
        #     t2_angle_dev = 0.1  # Angle deviation for P2 tangent
        #     t3_angle_dev = 0.1  # Angle deviation for P3 tangent

        p1_x_d = (0, 0)
        p1_y_d = (0, 0)
        p1_ang = (0, 0.1)

        p2_x_d = (1, 0.5)
        p2_y_d = (0.5, 0.1)
        p2_ang = (np.pi / 2, 0.1)

        p3_x_d = (0, 0.5)
        p3_y_d = (1, 0.1)
        p3_ang = (np.pi, 0.1)

        if self.point1 is None:
            self.point1 = Point(_pos(p1_x_d, p1_y_d), _rot([1, 0], _rand(*p1_ang)))
        if self.point2 is None:
            self.point2 = Point(_pos(p2_x_d, p2_y_d), _rot([1, 0], _rand(*p2_ang)))
        if self.point3 is None:
            self.point3 = Point(_pos(p3_x_d, p3_y_d), _rot([1, 0], _rand(*p3_ang)))

        self.base_points = PointGroup([self.point1, self.point2, self.point3])


    def _create_point_group(self):
        t = np.array([0, 1, 2])
        points   = np.array([self.point1.pos, self.point2.pos, self.point3.pos])
        tangents = np.array([self.point1.tangent, self.point2.tangent, self.point3.tangent])

        # Create Hermite spline for x and y
        spline_x = CubicHermiteSpline(t, points[:, 0], tangents[:, 0])
        spline_y = CubicHermiteSpline(t, points[:, 1], tangents[:, 1])

        s_vals = np.linspace(0, 2, self.num_points)

        self.point_group = PointGroup([Point(np.array([spline_x(s), spline_y(s)])) for s in s_vals])
        self.point_group.compute_curvature()


    def _create_fem_geometry(self):
        # convert point group to list of segments
        segments = []
        for i in range(1, len(self.point_group.points)):
            p_start = self.point_group.points[i - 1]
            p_end = self.point_group.points[i]
            segments.append(StraightSegment(p_start.pos, p_end.pos, n_subdivisions=1, name=f"Segment {i}"))

        p1 = self.base_points.points[-1].pos
        p2 = p1 + self.base_points.points[-1].tangent * 3
        p3 = p2 + self.base_points.points[-1].rotated(-np.pi/2).tangent  * 1
        p4 = p3 + self.base_points.points[-1].rotated( np.pi  ).tangent  * 6

        p0 = self.base_points.points[0].pos
        pm1 = p0 - self.base_points.points[0].tangent * 3
        pm2 = pm1 + self.base_points.points[0].rotated(-np.pi/2).tangent * 1
        pm3 = pm2 - self.base_points.points[0].rotated( np.pi  ).tangent  * 6

        segments.append(StraightSegment(p1, p2, n_subdivisions=20))
        segments.append(StraightSegment(p2, p3, n_subdivisions=10, name="LOAD"))
        segments.append(StraightSegment(p3, p4, n_subdivisions=20))
        segments.append(StraightSegment(p4, pm3, n_subdivisions=10))
        segments.append(StraightSegment(pm3, pm2, n_subdivisions=20))
        segments.append(StraightSegment(pm2, pm1, n_subdivisions=10, name="SUPP"))
        segments.append(StraightSegment(pm1, p0, n_subdivisions=20))

        geom_2d = Geometry.mesh_interior([SegmentGroup(segments)])
        geom_3d = geom_2d.as_second_order()
        geom_3d = geom_3d.extruded(1, -1)
        return geom_2d, geom_3d


    def _compute_ids(self):
        for i, p in enumerate(self.point_group.points):
            # go through all nodes in the fem geometry
            best_dist = np.inf
            best_id = None
            for node_id, node in enumerate(self.fem_geom_3d.nodes):
                if node is None:
                    continue
                dist = np.linalg.norm(np.asarray(node[:2]) - p.pos)
                if dist < best_dist:
                    best_dist = dist
                    best_id = node_id

            self.point_group.points[i].id = best_id

    def compute_stress(self, solver="./../../bin/FEMaster"):
        # compute fem geometry
        self.fem_geom_2d, self.fem_geom_3d = self._create_fem_geometry()
        self._compute_ids()

        # write input deck
        inp_deck = f"geometry{self.id}.inp"
        sol_deck = f"geometry{self.id}.res"
        self.fem_geom_3d.write_input_deck(inp_deck)
        with open(inp_deck, "a") as f:
            # material
            f.write("*MATERIAL, NAME=STEEL\n")
            f.write("*ELASTIC, TYPE=ISOTROPIC\n")
            f.write("200000., 0.3\n")
            f.write("*DENSITY\n")
            f.write("7.85e-9\n")
            f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL\n")
            # load
            f.write("*CLOAD, LOAD_COLLECTOR=LOAD\n")
            f.write("LOAD, 0., 1., 0.\n")
            # support
            f.write("*SUPPORT, SUPPORT_COLLECTOR=SUPP\n")
            f.write("SUPP, 0, 0, 0\n")
            # loadcase
            f.write("*LOAD CASE, TYPE=LINEAR STATIC\n")
            f.write("*LOAD\n")
            f.write("LOAD\n")
            f.write("*SUPPORT\n")
            f.write("SUPP\n")
            f.write("*END")

        # run solver and wait for completion
        os.system(f"{solver} {inp_deck}")

        sol = Solution.open(sol_deck)
        stress = sol.list_fields_reduced()['mises']()
        for i, p in enumerate(self.point_group.points):
            p.stress = stress[p.id]

        # remove all files
        os.remove(inp_deck)
        os.remove(sol_deck)

# generate 20 models in parallel
import multiprocessing

def generate_model(i):
    model = Model(id=i)
    model.compute_stress()
    return model

pool = multiprocessing.Pool(10)
models = pool.map(generate_model, range(200))

# sort models by max stress
models.sort(key=lambda m: max([p.stress for p in m.point_group.points]), reverse=True)

for m in models:
    print(f"Model {m.id} "
          f"- P1: {m.base_points.points[0].pos} "
          f"- P2: {m.base_points.points[1].pos} "
          f"- P3: {m.base_points.points[2].pos} "
          f"Max Stress: {max([p.stress for p in m.point_group.points]):.2f}")

# plot the top 5 models and worst models
fig, axs = plt.subplots(5, 2, figsize=(15, 15))
for i in range(5):
    for j in range(2):
        model = models[i + j * 195]
        points = np.array([p.pos for p in model.point_group.points])
        stress = [p.stress for p in model.point_group.points]

        axs[i, j].plot(points[:, 0], points[:, 1], 'o-')
        axs[i, j].set_title(f"Model {i + j * 195} - Max Stress: {max(stress):.2f}")
        axs[i, j].axis("equal")
        axs[i, j].grid(True)

plt.show()
