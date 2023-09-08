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
    for x in range(width + 1):
        for y in range(height + 1):
            for z in range(length + 1):
                nodes.append((x,y,z))
                node_lookup[(x, y, z)] = node_id
                node_id += 1

    elements = []
    element_id = 0
    for x in range(width):
        for y in range(height):
            for z in range(length):
                nodes_indices = [
                    node_lookup[(x  , y  , z  )], node_lookup[(x+1, y  , z  )],
                    node_lookup[(x+1, y+1, z  )], node_lookup[(x  , y+1, z  )],
                    node_lookup[(x  , y  , z+1)], node_lookup[(x+1, y  , z+1)],
                    node_lookup[(x+1, y+1, z+1)], node_lookup[(x  , y+1, z+1)]
                ]
                elements.append(nodes_indices)
                element_id += 1


    return nodes, elements, node_lookup

def get_special_node_ids(width, height, length, node_lookup):

    # Corner IDs
    corner1 = 0
    corner2 = width * (height + 1) * (length + 1)
    corner3 = height * (length + 1)
    corner4 = corner3 + corner2
    corner5 = length + corner1
    corner6 = length + corner2
    corner7 = length + corner3
    corner8 = length + corner4
    corners = [corner1, corner2, corner3, corner4, corner5, corner6, corner7, corner8]

    # Define edge midpoint coordinates
    edge_mids_coords = [
        ((width+1)//2, 0, 0),
        (width, (height+1)//2, 0),
        ((width+1)//2, height, 0),
        (0, (height+1)//2, 0),
        ((width+1)//2, 0, length),
        (width, (height+1)//2, length),
        ((width+1)//2, height, length),
        (0, (height+1)//2, length),
        (0, 0, (length+1)//2),
        (width, 0, (length+1)//2),
        (0, height, (length+1)//2),
        (width, height, (length+1)//2)
    ]

    # Get the IDs for edge midpoints
    edge_mids = [node_lookup[coords] for coords in edge_mids_coords]

    # Define face midpoint coordinates
    face_mids_coords = [
        ((width+1)//2, (height+1)//2, 0),
        ((width+1)//2, (height+1)//2, length),
        ((width+1)//2, 0, (length+1)//2),
        ((width+1)//2, height, (length+1)//2),
        (0, (height+1)//2, (length+1)//2),
        (width, (height+1)//2, (length+1)//2)
    ]

    # Get the IDs for face midpoints
    face_mids = [node_lookup[coords] for coords in face_mids_coords]

    # Define volume midpoint coordinate and get its ID
    volume_mid = node_lookup[((width+1)//2, (height+1)//2, (length+1)//2)]

    return corners, edge_mids, face_mids, volume_mid
