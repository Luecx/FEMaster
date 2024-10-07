from model import topo
from model import geometry
from model import tovtk

import numpy as np



################################################################################
# creating geometry
################################################################################
side_length = 10
side_divisions = 1

points = []
points.append([0, 0])
angle = 0
# example to optimise a hexagon
for i in range(4):
    for n in range(side_divisions):
        old_x = points[-1][0]
        old_y = points[-1][1]
        new_x = old_x + side_length / side_divisions * np.cos(angle)
        new_y = old_y + side_length / side_divisions * np.sin(angle)
        points.append([new_x, new_y])
    angle += np.pi / 2
points = np.asarray(points)



geom_2d = geometry.Geometry.mesh_interior(points, second_order=False, mesh_type=0)
geom_2d = geom_2d.extruded(10)
geom_2d.write_input_deck("hexagon_first_order.inp")

geom_2d_b = geometry.Geometry.mesh_interior(points, second_order=False,  mesh_type=2)
geom_2d_b = geom_2d_b.extruded(10)
geom_2d_b = geom_2d_b.as_second_order()
geom_2d_b.write_input_deck("hexagon_second_order.inp")

geom_2d_c = geometry.Geometry.mesh_interior(points, second_order=True,  mesh_type=2)
geom_2d_c.write_input_deck("hexagon_second_order_c1.inp")
geom_2d_c = geom_2d_c.extruded(1)
geom_2d_c.write_input_deck("hexagon_second_order_c2.inp")

tovtk.Converter("hexagon_second_order_c1.inp").convert()
tovtk.Converter("hexagon_second_order_c2.inp").convert()


# # opt = topo.Optimiser(input_deck="../../../../Desktop/stift_opt/01_input/designspace.inp",
# #                 desi_set="EALL",
# #                 loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                 output_folder="../../../../Desktop/stift_opt/opt6/")
# # opt = topo.Optimiser(input_deck="..\\..\\..\\..\\Desktop\\stift_opt\\01_input\\designspace.inp",
# #                      desi_set="EALL",
# #                      loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                      output_folder="..\\..\\..\\..\\Desktop\\stift_opt\\opt6\\")
# # opt = topo.Optimiser(input_deck="../../../../Desktop/whiteboard_opt/bauraum3.inp",
# #                      desi_set="DESI",
# #                      loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                      symmetries={'xy': 0.075},
# #                      output_folder="../../../../Desktop/whiteboard_opt/opt3/")
# # opt = topo.Optimiser(input_deck="..\\..\\..\\..\\Desktop\\whiteboard_opt\\bauraum3.inp",
# #                      desi_set="DESI",
# #                      loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                      output_folder="..\\..\\..\\..\\Desktop\\whiteboard_opt\\opt3\\")
# #
# # opt.set_solver(path="../bin/FEMaster.exe", method="direct", device="cpu", ncpus=6)
# # opt.set_exponent(3)
# # opt.set_target_density(0.1)
# # opt.set_filter(0.001)
# # opt.set_symmetry_radius(0.0005)
# # opt.set_move_limit(0.2)
# # # opt.start(50)
# # opt.load_it(18)
# # opt.plot(0.8)
# # opt.to_stl(0.8, "..\\..\\..\\..\\Desktop\\whiteboard_opt\\opt3\\result.stl")
#
# opt = topo.Optimiser(input_deck="res/haken_analysis_2/v2opt.inp",
#                     desi_set="EALL",
#                     loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                     output_folder="res/haken_analysis_2/opt2/")
# opt.set_solver(path="bin/FEMaster.exe", method="direct", device="cpu")
# opt.set_exponent(3)
# opt.set_target_density(0.4)
# opt.set_filter(3)
# opt.set_move_limit(0.2)
# opt.start(30)
#
# # opt = topo.Optimiser(input_deck="../res/topo_bridge_1/Job-2.inp",
# #                 desi_set="EALL",
# #                 loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                 output_folder="../res/topo_bridge_1/opt3/")
# # pt = topo.Optimiser(input_deck="..\\res\\topo_bridge_1\\Job-2.inp",
# #                     desi_set="EALL",
# #                     loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
# #                     output_folder="..\\res\\topo_bridge_1\\opt3\\")
# # opt.set_solver(path="../bin/FEMaster.exe", method="direct", device="cpu")
# # opt.set_exponent(3)
# # opt.set_target_density(0.2)
# # opt.set_filter(2)
# # opt.set_move_limit(0.2)
# # # opt.start(50)
# # opt.load_it(27)
# # opt.plot(0.3)
#
# #opt.load_it(50)
# #opt.plot(0.3)
#
# # from model import Solution
# #
# # sol = Solution.open("opt1/iterations/1/model.inp.res")
# # print(sol)
#
#
