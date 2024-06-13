import model.geometry as geometry
import model.topo as topo
import model.viewer as viewer
import model.solution as solution
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# generate a cuboid with dimensions LxWxH
def generate_cuboid(L, W, H, length_ratio=1):
    geom = geometry.Geometry()

    # lambda to index nodes using i, j, k
    node_id = lambda i, j, k: i * (W + 1) * (H + 1) + j * (H + 1) + k
    elem_id = lambda i, j, k: i * W * H + j * H + k

    # add nodes
    for i in range(L + 1):
        for j in range(W + 1):
            for k in range(H + 1):
                geom.add_node(node_id(i, j, k), i * length_ratio, j, k)

    # add elements
    for i in range(L):
        for j in range(W):
            for k in range(H):
                node_ids = [node_id(i, j, k),
                            node_id(i + 1, j, k),
                            node_id(i + 1, j + 1, k),
                            node_id(i, j + 1, k),
                            node_id(i, j, k + 1),
                            node_id(i + 1, j, k + 1),
                            node_id(i + 1, j + 1, k + 1),
                            node_id(i, j + 1, k + 1)]
                geom.add_element(elem_id(i, j, k), 'C3D8', node_ids)

    # add all the nodes which have x = 0 to a new set called "constraint"
    for j in range(W + 1):
        for k in range(H + 1):
            geom.add_node_to_set('CONSTRAINT_SET', node_id(0, j, k))

    # add the center node at the opposing end to a new set called "load"
    # if the amount of nodes in y and z is even, use the 4 center nodes
    # if the amount of nodes in y and z is odd, use the 1 center node
    if W % 2 == 1:
        geom.add_node_to_set('CENTER'   , node_id(L, W // 2, H // 2))
        geom.add_node_to_set('CENTER'   , node_id(L, W // 2 + 1, H // 2))
        geom.add_node_to_set('CENTER'   , node_id(L, W // 2, H // 2 + 1))
        geom.add_node_to_set('CENTER'   , node_id(L, W // 2 + 1, H // 2 + 1))

        geom.add_node_to_set('TOPCENTER', node_id(L, W, H // 2))
        geom.add_node_to_set('TOPCENTER', node_id(L, W, H // 2 + 1))
        geom.add_node_to_set('BOTCENTER', node_id(L, 0, H // 2))
        geom.add_node_to_set('BOTCENTER', node_id(L, 0, H // 2 + 1))
    else:
        geom.add_node_to_set('CENTER'   , node_id(L, W // 2, H // 2))

        geom.add_node_to_set('TOPCENTER', node_id(L, W , H // 2))
        geom.add_node_to_set('BOTCENTER', node_id(L, 0 , H // 2))
    geom.add_node_to_set('OUTER', node_id(L, 0, 0))
    geom.add_node_to_set('OUTER', node_id(L, W, 0))
    geom.add_node_to_set('OUTER', node_id(L, 0, H))
    geom.add_node_to_set('OUTER', node_id(L, W, H))
    # geom.add_node_to_set('LOAD_SET', node_id(L, 0, 0))
    # geom.add_node_to_set('LOAD_SET', node_id(L, W, H))
    # geom.add_node_to_set('LOAD_SET', node_id(L, 0, H))
    geom.change_to_second_order()

    return geom

def mean_density(x, L, W, H):
    temp = x.reshape((L, W, H))

    # Calculate the mean along the x-axis (axis 0)
    mean_temp = np.mean(temp[2 * L // 5 : 3 * L // 5], axis=0)

    # enforce symmetry along x dimension first, if its uneven, ignore the center
    mean_temp = (mean_temp + mean_temp[::-1,:]) / 2
    mean_temp = (mean_temp + mean_temp[:,::-1]) / 2

    # along the x-axis to restore the original shape.
    mean_expanded = np.repeat(mean_temp[np.newaxis, :, :], L, axis=0)

    return mean_expanded.flatten()

def create_file(N, length_ratio, moment, force, type='analysis'):
    W = N
    H = N
    L = 6 * N

    # compute moment force pair
    moment_force = moment / (H)

    name = "beam_{}_{}_{}_d".format(L, W, H)
    geom = generate_cuboid(L,W,H, length_ratio=length_ratio)
    geom.write_input_deck(f"{name}.inp")
    str = f"""
*MATERIAL, NAME=STEEL
*ELASTIC, TYPE=ISO
200000, 0.0
*DENSITY
7.85E-9
*SOLID SECTION, ELSET=EALL, MAT=STEEL

*CLOAD, LOAD_COLLECTOR=LOADS1
TOPCENTER, {moment_force}, {force / 2}, 0
BOTCENTER, {-moment_force}, {force / 2}, 0
*SUPPORT, SUPPORT_COLLECTOR=SUPPS
CONSTRAINT_SET, 0, 0, 0
"""

    if type == 'analysis':
        str += """
*LOADCASE, TYPE=LINEAR STATIC
*LOAD
LOADS1
*SUPPORT
SUPPS
*SOLVER, DEVICE=cpu, METHOD=indirect
*END
"""

    with open(f"{name}.inp", "a") as f:
        f.write(str)

    return f"{name}.inp", moment + force * L, force, (L, W, H)



for moment in [1, 100, 1000, 10000, 100000]:
    file_name, moment, force, (L, W, H) = create_file(14, 3, moment=moment, force=1, type='topo')
    opt = topo.Optimiser(input_deck=f"{file_name}",
                         solver_path="../bin/FEMaster.exe",
                         desi_set="EALL",
                         loadcases=[{"load_cols":["LOADS1"], "supp_cols": ["SUPPS"]}],
                         output_folder=f"./{file_name.split('.')[0]}_m{moment}_f{force}/",
                         filter_radius=0,
                         target_density=0.5,
                         move_limit=0.2,
                         min_density=0.001,
                         exponent=3,
                         method='indirect',
                         custom_density_adjustment=lambda x: mean_density(x, L, W, H))
    opt.start(50)
    # opt.load_it(10)
    # opt.plot(0.6)
    dens = opt.density[opt.desi_mask].reshape((L, W, H))[0,:,:]
    # write dens to a png file
    fig, ax = plt.subplots()

    plt.imshow(dens, cmap='gray', vmin=0, vmax=1)
    ax.axis('off')
    plt.savefig(f"./{moment}_density.png", bbox_inches='tight', pad_inches=0)
    plt.close(fig)