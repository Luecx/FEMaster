import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import Normalize

def plot_cuboid(model, densities, threshold=0.5, colormap='viridis'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Normalize the density values to be in the range [0, 1]
    norm = Normalize(vmin=min(densities), vmax=max(densities))

    for element, density in zip(model.elements, densities):
        if density >= threshold:
            # extract corner points of the current element
            corners = [model.nodes[i] for i in element]
            # define the 6 faces of the element
            faces = [
                [corners[j] for j in [0, 1, 2, 3]],
                [corners[j] for j in [4, 5, 6, 7]],
                [corners[j] for j in [0, 3, 7, 4]],
                [corners[j] for j in [1, 2, 6, 5]],
                [corners[j] for j in [0, 1, 5, 4]],
                [corners[j] for j in [2, 3, 7, 6]]
            ]

            # Use the colormap to determine the color based on density
            color = plt.cm.get_cmap(colormap)(norm(density))

            ax.add_collection3d(Poly3DCollection(faces, alpha=1, linewidths=1, edgecolors='b', facecolors=color))

    # Extract the minimum and maximum coordinates for each axis from nodes
    min_x, min_y, min_z = np.min(model.nodes, axis=0)
    max_x, max_y, max_z = np.max(model.nodes, axis=0)

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_zlim(min_z, max_z)

    # 8 corners of the cuboid
    corners = [(min_x, min_y, min_z), (max_x, min_y, min_z),
               (max_x, max_y, min_z), (min_x, max_y, min_z),
               (min_x, min_y, max_z), (max_x, min_y, max_z),
               (max_x, max_y, max_z), (min_x, max_y, max_z)]

    # 12 edges defined by pairs of corner indices
    edges = [(0, 1), (1, 2), (2, 3), (3, 0),
             (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]

    for edge in edges:
        p1, p2 = corners[edge[0]], corners[edge[1]]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='k')

    # B) Display nodes which are constrained:
    if model.supports:
        constrained_nodes = [model.nodes[id] for id in model.supports.keys()]
        for id in model.supports:
            point = model.nodes[id]
            type  = model.supports[id]
            tuple_to_str = lambda t: ''.join([char if val is not None else '_' for char, val in zip('xyz', t)])
            ax.text(point[0]+0.1, point[1], point[2], tuple_to_str(type), color='b')
        ax.scatter(*zip(*constrained_nodes), color='red', s=50, marker='o')  # Constrained nodes in red

    # C) Add a vector to the nodes that are loaded:
    if model.loads:
        for node_id, (dx, dy, dz) in model.loads.items():
            node = model.nodes[node_id]
            ax.quiver(node[0], node[1], node[2], dx, dy, dz, color='blue', length=0.5)  # Load vectors in blue



    axis_length = 3  # or another appropriate value
    origin = [-1, -1, -1]  # or another appropriate location

    # Draw the x, y, and z axis arrows
    ax.quiver(origin[0], origin[1], origin[2], axis_length, 0, 0, color='r', arrow_length_ratio=0.1)
    ax.quiver(origin[0], origin[1], origin[2], 0, axis_length, 0, color='g', arrow_length_ratio=0.1)
    ax.quiver(origin[0], origin[1], origin[2], 0, 0, axis_length, color='b', arrow_length_ratio=0.1)

    # Label the axes
    ax.text(origin[0] + axis_length + 0.1, origin[1], origin[2], "X", color='r')
    ax.text(origin[0], origin[1] + axis_length + 0.1, origin[2], "Y", color='g')
    ax.text(origin[0], origin[1], origin[2] + axis_length + 0.1, "Z", color='b')


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

# Example usage:
# Replace nodes, elements, and densities with your actual data
# plot_cuboid(nodes, elements, densities)
