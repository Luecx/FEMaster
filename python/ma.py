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
L = 100
W = 21
H = 3


def mean_density(x):
    temp = x.reshape((L, W, H))

    # Calculate the mean along the x-axis (axis 0)
    mean_temp = np.mean(temp[2 * L // 5 : 4 * L // 5], axis=0)

    # enforce symmetry along x dimension first, if its uneven, ignore the center
    mean_temp = (mean_temp + mean_temp[::-1,:]) / 2
    mean_temp = (mean_temp + mean_temp[:,::-1]) / 2


# The output should have the same shape as the input, so we need to repeat the mean values
    # along the x-axis to restore the original shape.
    mean_expanded = np.repeat(mean_temp[np.newaxis, :, :], L, axis=0)

    return mean_expanded.flatten()

def run_target(target, length_ratio):

    name = "beam_{}_{}_{}_d".format(L, W, H)

    # generate a cuboid with dimensions LxWxH
    geom = generate_cuboid(L,W,H, length_ratio=length_ratio)
    geom.write_input_deck(f"{name}.inp")

    # append the loads to the input deck:

    str = """
    *MATERIAL, NAME=STEEL
    *ELASTIC, TYPE=ISO
    200000, 0.0
    *DENSITY
    7.85E-9
    *SOLID SECTION, ELSET=EALL, MAT=STEEL
    
    *CLOAD, LOAD_COLLECTOR=LOADS1
    TOPCENTER, 0, 1, 0
    BOTCENTER, 0, 1, 0
    ***CLOAD, LOAD_COLLECTOR=LOADS2
    **LOAD_SET, 0, 0, 1
    *SUPPORT, SUPPORT_COLLECTOR=SUPPS
    CONSTRAINT_SET, 0, 0, 0
    """

    with open(f"{name}.inp", "a") as f:
        f.write(str)

    opt = topo.Optimiser(input_deck=f"{name}.inp",
                    solver_path="../bin/FEMaster.exe",
                    desi_set="EALL",
                         loadcases=[{"load_cols":["LOADS1"], "supp_cols": ["SUPPS"]}],
                    output_folder=f"./{name}/",
                    filter_radius=0,
                    target_density=target,
                         move_limit=0.2,
                     min_density=0.0001,
                         exponent=3,
                    method='direct',
                     custom_density_adjustment=mean_density)

    opt.start(1)
    return None

    opt.load(f"./{name}/iterations/30/model.dat")
    sol = solution.Solution.open(f"./{name}/iterations/30/model.inp.res")

    # elem_mask = np.zeros(L * W * H +1 , dtype=bool)
    # for i in range(L // 2):
    #     for j in range(W):
    #         for k in range(H):
    #             elem_mask[i * W * H + j * H + k + 1] = True
    #
    # stress = sol.loadcases['1']['STRESS']
    # strain = sol.loadcases['1']['STRAIN']
    # displc = sol.list_fields()['displacement']()
    # dispxy = sol.list_fields()['displacement_xyz']()
    #
    #
    # # element wise product and sum over
    # # the last axis (the stress components)
    # comp = np.sum(stress * strain, axis=1)
    # print(comp.shape)
    # print(np.min(comp), np.max(comp))
    #
    # dens_mask = opt.density > 0.1
    #
    # node_cords = []
    # for node in opt.geometry.nodes:
    #     if node is None:
    #         continue
    #     node_cords.append(list(node))
    # node_cords = np.asfarray(node_cords)
    #
    # max_disp_val = np.max(np.abs(dispxy))
    # max_cord_dist = np.max(np.max(node_cords, axis=1) - np.min(node_cords, axis=1))
    #
    # v = viewer.Viewer()
    # v.set_geometry(opt.geometry)
    # v.set_element_mask(elem_mask & dens_mask)
    # v.set_data(type='node', data=comp)
    # #v.set_data(type='element', data=opt.density)
    # #v.set_data(type='element', data=sol.list_fields()['compliance_raw']())
    # v.set_displacement(dispxy * (max_cord_dist / max_disp_val) / 10)
    # v.set_data_range(0,1,True)
    # # v.set_data_range(-0.0001, 0.0001, False)
    # v.set_colorscheme('jet')
    # v.coordinate_system()
    # v.set_boundaries()
    # v.set_grid_xy()
    # v.mainloop()

    def animate(i):
        # Load the data for each iteration
        opt.load(f"./{name}/iterations/{i}/model.dat")

        # sol = solution.Solution.open(f"./{name}/iterations/{i}/model.inp.res")
        #grd = sol.list_fields()['compliance_raw']()[opt.desi_mask]

        #opt.load(f"./{name}/iterations/{i}/model.dat")
        #dens = opt.density[opt.desi_mask].reshape((L, W, H))[0,:,:]

        grd = opt.density[opt.desi_mask]

        # normalise grd
        grd = (grd - np.min(grd)) / (np.max(grd) - np.min(grd))

        grd = grd.reshape((L, W, H))[L // 2,:,:]

        # Clear the current image and plot the new one
        ax.clear()
        im = ax.imshow(grd, cmap='gray', vmin=0, vmax=1)
        return im,

    # Create the animation
    fig, ax = plt.subplots()
    ax.set_xlim(0, W)  # Set the x-axis limits
    ax.set_ylim(0, L)  # Set the y-axis limits
    ani = FuncAnimation(fig, animate, frames=range(1, 40), interval=200, blit=True)
    plt.show()

    return opt.density[opt.desi_mask].reshape((L, W, H))[0,:,:]


run_target(0.3, 100)
#density = run_target(0.5, 100)

# for length_ratio in [3, 10, 100, 1000]:
#     # run for density increments of 0.05 and save the results as images using matplotlib
#     for target in [0.1, 0.3, 0.5]:
#         density = run_target(target, length_ratio)
#
#         # Create a figure and an axes object
#         fig, ax = plt.subplots()
#
#         # Display the image in grayscale, and remove axes and titles
#         ax.imshow(density, cmap='gray', vmin=0, vmax=1)
#         ax.axis('off')  # This removes the axis
#
#         # Save the figure with a tight layout and no padding
#         plt.savefig(f"den_{target}_len_{length_ratio}.png", bbox_inches='tight', pad_inches=0)
#
#         # Close the figure to free memory
#         plt.close(fig)

# opt.load(f"./{name}/iterations/40/model.dat")

# sol = solution.Solution.open(f"./{name}/iterations/40/model.inp.res")
#
# v = viewer.Viewer()
# v.set_geometry(geometry.Geometry.read_input_deck(f"./{name}/iterations/40/model.inp"))
# v.set_element_mask(opt.density > 0.1)
# v.set_data(type='element', data=opt.density)
#
# print(sol.list_fields()['mises']())
# #v.set_data(sol.list_fields()['mises']())
# v.set_colorscheme('jet')
# v.set_data_range(0,1, 1)
# v.set_displacement(sol.list_fields()['displacement_xyz']() * 1000)
#
# v.coordinate_system()
# v.set_grid_xy()
# v.mainloop()


# opt.start(100)
#opt.load("./../beam_10x10/iterations/100/model.dat")
#opt.plot(threshold=0.3)


