import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_cuboid(nodes, elements, densities, threshold=0.5):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for element, density in zip(elements, densities):
        if density >= threshold:
            # extract corner points of the current element
            corners = [nodes[i] for i in element]
            # define the 6 faces of the element
            faces = [
                [corners[j] for j in [0, 1, 2, 3]],
                [corners[j] for j in [4, 5, 6, 7]],
                [corners[j] for j in [0, 3, 7, 4]],
                [corners[j] for j in [1, 2, 6, 5]],
                [corners[j] for j in [0, 1, 5, 4]],
                [corners[j] for j in [2, 3, 7, 6]]
            ]
            ax.add_collection3d(Poly3DCollection(faces, alpha=1, linewidths=1, edgecolors='b'))

    # Extract the minimum and maximum coordinates for each axis from nodes
    min_x = min(node[0] for node in nodes)
    max_x = max(node[0] for node in nodes)
    min_y = min(node[1] for node in nodes)
    max_y = max(node[1] for node in nodes)
    min_z = min(node[2] for node in nodes)
    max_z = max(node[2] for node in nodes)

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_zlim(min_z, max_z)

    # Setting the aspect ratio
    range_x = max_x - min_x
    range_y = max_y - min_y
    range_z = max_z - min_z
    ax.set_box_aspect([range_x, range_y, range_z])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.axis('off')  # Turn off the axis
    plt.show()