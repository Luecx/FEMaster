import numpy as np

def compute_element_midpoints(nodes, elements):
    midpoints = []
    for node_ids in elements:
        coords = [nodes[node_id] for node_id in node_ids]
        midpoint = np.mean(coords, axis=0)
        midpoints.append(tuple(midpoint))
    return midpoints

def generate_geometry(width, height, length):

    nodes = []
    node_lookup = {}  # This is our temporary lookup
    node_id = 0
    for z in range(length + 1):
        for x in range(width + 1):
            for y in range(height + 1):
                nodes.append((x,y,z))
                node_lookup[(x, y, z)] = node_id
                node_id += 1

    elements = []
    element_id = 0
    for z in range(length):
        for x in range(width):
            for y in range(height):
                nodes_indices = [
                    node_lookup[(x  , y  , z  )], node_lookup[(x+1, y  , z  )],
                    node_lookup[(x+1, y+1, z  )], node_lookup[(x  , y+1, z  )],
                    node_lookup[(x  , y  , z+1)], node_lookup[(x+1, y  , z+1)],
                    node_lookup[(x+1, y+1, z+1)], node_lookup[(x  , y+1, z+1)]
                ]
                elements.append(nodes_indices)
                element_id += 1


    return nodes, elements


def get_special_node_ids(width, height, length):

    # Corner IDs
    corner1 = 0
    corner2 = width
    corner3 = width * (height + 1)
    corner4 = corner3 + width

    # For z=length
    offset_z = (width + 1) * (height + 1) * length
    corner5 = offset_z + corner1
    corner6 = offset_z + corner2
    corner7 = offset_z + corner3
    corner8 = offset_z + corner4
    corners = [corner1, corner2, corner3, corner4, corner5, corner6, corner7, corner8]

    # Mid-surface IDs
    mid_surfaces = [
        (width // 2),
        (width + 1) * (height // 2),
        (width + 1) * (height + 1) * (length // 2),
        (width + 1) * (height + 1) * (length // 2) + (width // 2),
        (width + 1) * (height + 1) * (length // 2) + (width + 1) * (height // 2),
        (width + 1) * (height + 1) * (length // 2) + width
    ]

    return corners, mid_surfaces