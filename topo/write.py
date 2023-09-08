import model
from stl import mesh
import numpy as np
import os
from plyfile import PlyData, PlyElement
import matplotlib.cm as cm
def mesh_to_stl(nodes, elements, densities, filename, threshold=0.5):
    # List to hold the triangles
    all_triangles = []

    for element, density in zip(elements, densities):
        if density >= threshold:
            corners = [nodes[i] for i in element]

            # Define the two triangles for each face of the cuboid
            triangles_for_faces = [
                [corners[j] for j in [0, 1, 2]], [corners[j] for j in [0, 2, 3]],
                [corners[j] for j in [4, 5, 6]], [corners[j] for j in [4, 6, 7]],
                [corners[j] for j in [0, 3, 7]], [corners[j] for j in [0, 7, 4]],
                [corners[j] for j in [1, 2, 6]], [corners[j] for j in [1, 6, 5]],
                [corners[j] for j in [0, 1, 5]], [corners[j] for j in [0, 5, 4]],
                [corners[j] for j in [2, 3, 7]], [corners[j] for j in [2, 7, 6]]
            ]

            all_triangles.extend(triangles_for_faces)

    # Create directories if they don't exist
    directory = os.path.dirname(filename)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

    # Create a new mesh object with the triangles
    mesh_object = mesh.Mesh(np.array(all_triangles, dtype=mesh.Mesh.dtype))
    for i, triangle in enumerate(all_triangles):
        mesh_object.vectors[i] = triangle

    # Write the mesh to an STL file
    mesh_object.save(filename)




def mesh_to_ply(nodes, elements, densities, filename, threshold=0.5):
    # For colormap interpolation
    min_density = min(densities)
    max_density = max(densities)
    norm = cm.colors.Normalize(vmin=min_density, vmax=max_density, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)  # using the viridis colormap

    # Prepare lists to hold vertex and face data
    all_vertices = []
    all_faces = []

    vertex_index = 0  # This keeps track of the current vertex index
    for element, density in zip(elements, densities):
        if density >= threshold:
            corners = [nodes[i] for i in element]
            color = (np.array(mapper.to_rgba(density)[:3]) * 255).astype(np.uint8)

            # For each cube, add the 8 vertices
            for corner in corners:
                all_vertices.append(tuple(corner) + tuple(color))

            # Define the two triangles for each face of the cuboid
            face_indices = [
                [3, 2, 1], [3, 1, 0],  # bottom
                [6, 7, 4], [6, 4, 5],  # top
                [3, 0, 4], [3, 4, 7],  # front
                [0, 1, 5], [0, 5, 4],  # right
                [1, 2, 6], [1, 6, 5],  # back
                [2, 3, 7], [2, 7, 6]   # left
            ]

            for face in face_indices:
                all_faces.append([vertex_index + idx for idx in face])

            vertex_index += 8

    # Create the PlyElements
    vertex_element = PlyElement.describe(
        np.array(all_vertices, dtype=[
            ('x', 'f4'), ('y', 'f4'), ('z', 'f4'),
            ('red', 'u1'), ('green', 'u1'), ('blue', 'u1')
        ]),
        'vertex'
    )

    # Create a structured numpy array with a separate field for vertex indices.
    faces_with_count = np.zeros(len(all_faces), dtype=[('vertex_indices', 'i4', (3,))])

    for idx, face in enumerate(all_faces):
        faces_with_count[idx]['vertex_indices'] = face

    face_element = PlyElement.describe(faces_with_count, 'face')


# Create directories if they don't exist
    directory = os.path.dirname(filename)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

    # Write to the ply file
    ply_data = PlyData([vertex_element, face_element], text=True)
    ply_data.write(filename)

def generate_model(model, densities):
    output = []
    output.append("*NODE")
    for node_id, (x, y, z) in enumerate(model.nodes):
        output.append(f"{node_id}, {x}, {y}, {z}")  # Adding 1 because node_id starts from 0

    output.append("\n*ELEMENT, TYPE=C3D8")
    for elem_id, node_ids in enumerate(model.elements):
        output.append(f"{elem_id}, " + ", ".join(map(str, [nid for nid in node_ids])))  # Adding 1 because indexing starts from 0

    output.append(f"\n*MATERIAL, NAME=MAT1")
    output.append(f"*ELASTIC, TYPE=ISO")
    output.append(f"{model.youngs},{model.nu}")
    output.append(f"*SOLID SECTION, ELSET=EALL, MAT=MAT1")

    output.append("\n*SUPPORT, SUPPORT_COLLECTOR=END_SUPPORTS")
    for node_id, values in model.supports.items():
        output.append(f"{node_id}"
                      f", {' ' if model.supports[node_id][0] is None else model.supports[node_id][0]}"
                      f", {' ' if model.supports[node_id][1] is None else model.supports[node_id][1]}"
                      f", {' ' if model.supports[node_id][2] is None else model.supports[node_id][2]}")

    output.append("\n*CLOAD, LOAD_COLLECTOR=END_LOADS")
    for node_id, values in model.loads.items():
        output.append(f"{node_id}"
                      f", {' ' if model.loads[node_id][0] is None else model.loads[node_id][0]}"
                      f", {' ' if model.loads[node_id][1] is None else model.loads[node_id][1]}"
                      f", {' ' if model.loads[node_id][2] is None else model.loads[node_id][2]}")

    if densities is not None:
        output.append(f"\n*LOADCASE, TYPE= LINEAR STATIC TOPO")
    else:
        output.append(f"\n*LOADCASE, TYPE= LINEAR STATIC")
    output.append(f"*SUPPORT \nEND_SUPPORTS")
    output.append(f"*LOAD \nEND_LOADS")

    if densities is not None:
        output.append("\n*DENSITY")
        for i, density in enumerate(densities):
            output.append(f"{i},{density}")

    output.append(f"\n*EXPONENT\n{model.exponent}")
    output.append(f"\n*SOLVER, METHOD={'DIRECT' if model.direct else 'INDIRECT'}, DEVICE={'CPU' if model.cpu else 'GPU'}\n*END")

    return '\n'.join(output)
