import numpy as np
import model
from scipy.optimize import minimize
import time
def compute_threshold_bisection(densities, target_density, tol=1e-6):
    """
    Compute the threshold above which the number of elements displayed would equal the target density using bisection.

    Parameters:
    - densities: array of element densities
    - target_density: desired target density
    - tol: tolerance for convergence (default is 1e-6)

    Returns:
    - threshold: the computed threshold
    """
    # Define lower and upper bounds for the threshold
    lower_bound, upper_bound = 0, 1

    # Target number of elements to display
    target_num_elements = target_density * len(densities)

    # Perform the bisection loop
    while upper_bound - lower_bound > tol:
        # Guess a threshold in between
        guess_threshold = (lower_bound + upper_bound) / 2

        # Count the number of elements above this threshold
        count = np.sum(densities > guess_threshold)

        # Adjust bounds based on the count
        if count > target_num_elements:
            lower_bound = guess_threshold
        else:
            upper_bound = guess_threshold

    return (lower_bound + upper_bound) / 2

model = model.Model(width=10, height=10, length=100)
corners, edge_mids, face_mids, volume_mid = model.get_special_node_ids()

model.set_support(corners[0], 0, 0, 0)
model.set_support(corners[1], 0, 0, 0)
model.set_support(corners[2], 0, 0, 0)
model.set_support(corners[3], 0, 0, 0)

model.set_force(face_mids[1], 0, 1000, 0)
# model.set_force(edge_mids[6], 0, 5, 0)
# model.set_force(corners[5], 1, 0, 0)
# model.set_force(corners[4], 0, -1, 0)
# model.set_force(corners[6], -1, 0, 0)
# model.set_force(corners[7], 0, 1, 0)
# model.set_force(corners[5], 0, 5, 0)
# model.set_force(corners[4], 0, 5, 0)
# model.set_force(corners[6], 0, 5, 0)
# model.set_force(corners[7], 0, 5, 0)

model.set_material(youngs=210000, nu=0.0)

model.add_symmetry("xz")
model.add_symmetry("yz")
model.add_symmetry("z") # rotational symmetry

model.set_exponent(2.5)
model.set_proximity_radius(3)
model.set_solver(use_cpu=False, use_direct=False)

model.write(None, "beam_mini.inp")

target_density = 1
n_elements = model.width * model.height * model.length

# x = np.ones(n_elements) * target_density
# loop = 0
# change = 1
# move = 0.2
#
# while change > 0.01 and loop < 50:
#     loop += 1
#
#     start_time = time.time()
#     opt_res = model.run(densities=x)
#     mid_time = time.time()
#     grad = opt_res["DENS_GRAD"]
#     comp = opt_res["COMPLIANCE_ADJ"]
#
#     x = np.array(x)
#     grad = np.array(grad)
#     grad = model.filter(grad)
#     grad = model.enforce_symmetry(grad)
#
#     # OC update of densities
#     l1, l2 = 1e-30, 1e30
#     xnew = np.zeros_like(x)
#     while (l2 - l1) / l1 > 1e-6:
#         lmid = 0.5 * (l2 + l1)
#         xnew[:] = np.maximum(0.01, np.maximum(x - move,  np.minimum(1, np.minimum(x + move, x*np.sqrt(-grad/(lmid))))))
#         xnew = model.enforce_symmetry(xnew)
#
#         if np.mean(xnew) - target_density > 0:
#             l1 = lmid
#         else:
#             l2 = lmid
#     change = np.linalg.norm(x - xnew)
#     x[:] = xnew
#
#     threshold = compute_threshold_bisection(x, target_density)
#
#     model.write_stl(f"{model.path}/stl/res_{loop}.stl", densities=x, threshold=threshold)
#     model.write_ply(f"{model.path}/ply/res_{loop}.ply", densities=x, threshold=threshold)
#     end_time = time.time()
#     print(f"Iteration {loop}: \n\t"
#           f"Objective  : {np.sum(comp):.5g}\n\t"
#           f"Change     : {change:.5g}\n\t"
#           f"Constraint : {np.mean(xnew):.5g}\n\t"
#           f"Threshold  : {threshold:.5g}\n\t"
#           f"Solver time: {mid_time - start_time:.5g}s\n\t"
#           f"Self time  : {end_time - mid_time:.5g}s")
#
# x = model.enforce_symmetry(x)
#
# threshold = compute_threshold_bisection(x, target_density)
#
# model.plot_cuboid(densities=x, threshold=threshold)