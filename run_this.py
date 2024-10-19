from geometry import *
from topopt import *

import numpy as np
inner_radius = 5
outer_radius = 15
thickness = 12
distance = 100


# computing the angle at which the outer circle intersects the strip which has a thickness
angle = np.arcsin(thickness / 2 / outer_radius)

offset_x = np.cos(angle) * outer_radius
offset_y = thickness / 2

outer = SegmentGroup([
    StraightSegment([offset_x, offset_y], [distance - offset_x, offset_y], subdivisions=100, name="TOP"),
    CircleSegment  (start_point=[distance - offset_x, offset_y], center=[distance,0], angles=-2*np.pi+2*angle, subdivisions=30),
    StraightSegment([distance - offset_x, -offset_y], [offset_x, -offset_y], subdivisions=100),
    CircleSegment  (start_point=[offset_x, -offset_y], center=[0,0], angles=-2*np.pi+2*angle, subdivisions=30),
])

inner1 = SegmentGroup([
    CircleSegment  (start_point=[inner_radius, 0], center=[0,0], angles=2*np.pi, subdivisions=30, name="INNER_1"),
])

inner2 = SegmentGroup([
    CircleSegment  (start_point=[distance + inner_radius, 0], center=[distance,0], angles=2*np.pi, subdivisions=30, name="INNER_2"),
])

geom = Geometry.mesh_interior([outer, inner1, inner2])
geom = geom.extruded(1, spacing=0.1)
geom = geom.as_second_order()

# create two new center nodes
id1 = geom.add_node(x=0, y=0, z=0)
id2 = geom.add_node(x=distance, y=0, z=0)

# create new sets and call them master
geom.add_node_to_set("MASTER_1", id1)
geom.add_node_to_set("MASTER_2", id2)

geom.write_input_deck("my_deck.inp")

# open file and append material and stuff
with open("my_deck.inp", 'a') as myfile:
    myfile.write("*Material, name=mat\n")
    myfile.write("*Elastic, type=iso\n")
    myfile.write("210000,0.3\n")
    myfile.write("*Density\n")
    myfile.write("7.85e-9\n")
    myfile.write("\n")
    myfile.write("*Solid Section, elset=EALL, material=mat\n")
    myfile.write("\n")
    myfile.write("*Coupling, type=kinematic, master=MASTER_1, slave=INNER_1\n")
    myfile.write("1,1,1,1,1,1\n")
    myfile.write("*Coupling, type=kinematic, master=MASTER_2, slave=INNER_2\n")
    myfile.write("1,1,1,1,1,1\n")
    myfile.write("\n")
    myfile.write("*Cload, load_collector=loads\n")
    myfile.write("TOP, 0, -1, 0\n")
    myfile.write("\n")
    myfile.write("*Support, support_collector=supp\n")
    myfile.write("INNER_1, 0, 0, 0,  , 0, 0\n")
    myfile.write("INNER_2,  , 0, 0,  , 0, 0\n")
    myfile.write("\n")

# create optimization problem
optimiser = Optimiser(
    input_deck="my_deck.inp",
    desi_set="EALL",
    loadcases=[{'load_cols': ['loads'], 'supp_cols': [1.0]}],
    output_folder="./my_deck_opt/",
    solver_path="../bin/FEMaster.exe",
    method='direct',
    device='cpu',
    exponent=2.5,
    min_density=0.01,
    target_density=0.2,
    filter_radius=3)

optimiser.start(50)


# outer = SegmentGroup([StraightSegment([0, 0], [5, 0], subdivisions=40),
#                       StraightSegment([5, 0], [5, 5], subdivisions=40),
#                       StraightSegment([5, 5], [0, 5], subdivisions=40),
#                       StraightSegment([0, 5], [0, 0], subdivisions=40)])
#
# inner = SegmentGroup([StraightSegment([2, 2], [3, 2], subdivisions=10),
#                       StraightSegment([3, 2], [3, 3], subdivisions=10),
#                       CurvedSegment  ([3, 3], [2, 3], mid=[2.5, 3.5], subdivisions=10),
#                       StraightSegment([2, 3], [2, 2], subdivisions=10, name="MY_CUSTOM_NAME_1")])
#
# # from 4 to 4.5
# inner2 = SegmentGroup([StraightSegment([4, 4], [4.5, 4], subdivisions=10, name="HERE_I_WANT_LOADS"),
#                        StraightSegment([4.5, 4], [4.5, 4.5], subdivisions=10),
#                        StraightSegment([4.5, 4.5], [4, 4.5], subdivisions=10),
#                        StraightSegment([4, 4.5], [4, 4], subdivisions=10)])
#
# geom = Geometry.mesh_interior([outer, inner, inner2])
# geom = geom.extruded(5, spacing=0.1)
# geom = geom.as_second_order()
# geom.write_input_deck("my_deck.inp")
