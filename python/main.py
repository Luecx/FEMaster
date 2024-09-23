
from model import topo


# opt = topo.Optimiser(input_deck="../../../../Desktop/stift_opt/01_input/designspace.inp",
#                 desi_set="EALL",
#                 loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                 output_folder="../../../../Desktop/stift_opt/opt6/")
# opt = topo.Optimiser(input_deck="..\\..\\..\\..\\Desktop\\stift_opt\\01_input\\designspace.inp",
#                      desi_set="EALL",
#                      loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                      output_folder="..\\..\\..\\..\\Desktop\\stift_opt\\opt6\\")
# opt = topo.Optimiser(input_deck="../../../../Desktop/whiteboard_opt/bauraum3.inp",
#                      desi_set="DESI",
#                      loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                      symmetries={'xy': 0.075},
#                      output_folder="../../../../Desktop/whiteboard_opt/opt3/")
opt = topo.Optimiser(input_deck="..\\..\\..\\..\\Desktop\\whiteboard_opt\\bauraum3.inp",
                     desi_set="DESI",
                     loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
                     output_folder="..\\..\\..\\..\\Desktop\\whiteboard_opt\\opt3\\")

opt.set_solver(path="../bin/FEMaster.exe", method="direct", device="cpu", ncpus=6)
opt.set_exponent(3)
opt.set_target_density(0.1)
opt.set_filter(0.001)
opt.set_symmetry_radius(0.0005)
opt.set_move_limit(0.2)
# opt.start(50)
opt.load_it(18)
opt.plot(0.8)
opt.to_stl(0.8, "..\\..\\..\\..\\Desktop\\whiteboard_opt\\opt3\\result.stl")

# opt = topo.Optimiser(input_deck="../res/topo_bridge_1/Job-2.inp",
#                 desi_set="EALL",
#                 loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                 output_folder="../res/topo_bridge_1/opt3/")
# pt = topo.Optimiser(input_deck="..\\res\\topo_bridge_1\\Job-2.inp",
#                     desi_set="EALL",
#                     loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
#                     output_folder="..\\res\\topo_bridge_1\\opt3\\")
# opt.set_solver(path="../bin/FEMaster.exe", method="direct", device="cpu")
# opt.set_exponent(3)
# opt.set_target_density(0.2)
# opt.set_filter(2)
# opt.set_move_limit(0.2)
# # opt.start(50)
# opt.load_it(27)
# opt.plot(0.3)

#opt.load_it(50)
#opt.plot(0.3)

# from model import Solution
#
# sol = Solution.open("opt1/iterations/1/model.inp.res")
# print(sol)


