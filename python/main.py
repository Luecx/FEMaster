from model import topo

opt = topo.Optimiser(input_deck="../res/topo_bridge_1/Job-2.inp",
                desi_set="EALL",
                loadcases=[{"load_cols": ["LOADS"], "supp_cols": ["SUPPS"]}],
                output_folder="../res/topo_bridge_1/opt/",
                     symmetries={"xy": 10}, symmetry_radius=1)
opt.set_solver(path="../bin/FEMaster.exe", method="direct", device="cpu")
opt.set_exponent(3)
opt.set_target_density(0.2)
opt.set_filter(3)
opt.set_move_limit(0.2)
opt.start(iterations=50)
# opt.load_it(6)
# opt.plot(0.3)

# from model import Solution
#
# sol = Solution.open("opt1/iterations/1/model.inp.res")
# print(sol)