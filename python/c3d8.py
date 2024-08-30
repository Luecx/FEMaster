import numpy as np
import sympy as sp

# Nodes and elements definition
nodes = {
    1: (-1, -1, -1),
    2: ( 1, -1, -1),
    3: ( 1,  1, -1),
    4: (-1,  1, -1),
    5: (-1, -1,  1),
    6: ( 1, -1,  1),
    7: ( 1,  1,  1),
    8: (-1,  1,  1),
}

elem_ids = [1, 2, 3, 4, 5, 6, 7, 8]

# Normal and reduced integration points
normal_integration = [
    ( 1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), 1.0),
    ( 1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), 1.0),
    ( 1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), 1.0),
    ( 1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), 1.0),
    (-1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), 1.0),
    (-1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), 1.0),
    (-1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0),  1.0/np.sqrt(3.0), 1.0),
    (-1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), -1.0/np.sqrt(3.0), 1.0)
]

reduced_integration = [(0, 0, 0, 8)]

# Material properties
material = {
    'E': 1,
    'nu': 0.0
}
def material_matrix_deviatoric_volumetric(material):
    E = material['E']
    nu = material['nu']

    # Lame parameters
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    # Total material matrix
    D_total = np.array(
        [[E, 0, 0, 0, 0, 0],
         [0, E, 0, 0, 0, 0],
         [0, 0, E, 0, 0, 0],
         [0, 0, 0, mu, 0, 0],
         [0, 0, 0, 0, mu, 0],
         [0, 0, 0, 0, 0, mu]])

    # Volumetric part of the material matrix
    K = lambda_ + (2 / 3) * mu
    D_vol = np.array(
        [[K, lambda_, lambda_, 0, 0, 0],
         [lambda_, K, lambda_, 0, 0, 0],
         [lambda_, lambda_, K, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]])

    # Deviatoric part of the material matrix
    D_dev = D_total - D_vol

    return D_dev, D_vol


def shape_functions():
    r, s, t = sp.symbols('x y z')
    N1 = (1 / 8) * (1 - r) * (1 - s) * (1 - t)
    N2 = (1 / 8) * (1 + r) * (1 - s) * (1 - t)
    N3 = (1 / 8) * (1 + r) * (1 + s) * (1 - t)
    N4 = (1 / 8) * (1 - r) * (1 + s) * (1 - t)
    N5 = (1 / 8) * (1 - r) * (1 - s) * (1 + t)
    N6 = (1 / 8) * (1 + r) * (1 - s) * (1 + t)
    N7 = (1 / 8) * (1 + r) * (1 + s) * (1 + t)
    N8 = (1 / 8) * (1 - r) * (1 + s) * (1 + t)
    shape_funcs = [N1, N2, N3, N4, N5, N6, N7, N8]
    return shape_funcs

def stress_strain_b(shape_functions_derivatives, integration_point):
    x, y, z, _ = integration_point
    B = np.zeros((6, 24))
    dN_subs = shape_functions_derivatives.subs({'x': x, 'y': y, 'z': z})
    for i in range(8):
        B[0, 3 * i] = dN_subs[i, 0]
        B[1, 3 * i + 1] = dN_subs[i, 1]
        B[2, 3 * i + 2] = dN_subs[i, 2]
        B[3, 3 * i] = dN_subs[i, 1] / 2
        B[3, 3 * i + 1] = dN_subs[i, 0] / 2
        B[4, 3 * i + 1] = dN_subs[i, 2] / 2
        B[4, 3 * i + 2] = dN_subs[i, 1] / 2
        B[5, 3 * i] = dN_subs[i, 2] / 2
        B[5, 3 * i + 2] = dN_subs[i, 0] / 2
    return B

def stiffess_at_integration_point(shape_functions_derivatives, integration_point, D):
    B = stress_strain_b(shape_functions_derivatives, integration_point)
    return np.dot(np.dot(B.T, D), B)

def shape_function_derivatives(shape_functions):
    derivatives_matrix = []
    for N in shape_functions:
        derivatives = [sp.diff(N, var) for var in ('x', 'y', 'z')]
        derivatives_matrix.append(derivatives)
    return sp.Matrix(derivatives_matrix)

def stiffness(shape_function_derivatives):
    D_dev, D_vol = material_matrix_deviatoric_volumetric(material)
    total_stiffness = np.zeros((24, 24))
    for integration_point in normal_integration:
        total_stiffness += integration_point[3] * stiffess_at_integration_point(shape_function_derivatives, integration_point, D_dev)
    for integration_point in reduced_integration:
        total_stiffness += integration_point[3] * stiffess_at_integration_point(shape_function_derivatives, integration_point, D_vol)
    return total_stiffness

shape_functions = shape_functions()
shape_functions_derivatives = shape_function_derivatives(shape_functions)

constraints = np.array([True, True, True,
                        False, False, False,
                        False, False, False,
                        True, True, True,
                        True, True, True,
                        False, False, False,
                        False, False, False,
                        True, True, True,])

loads = np.array([0, 0, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 0, 0,
                  0, 0, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 0, 0,])

np.set_printoptions(precision=5, suppress=True, linewidth=20000)

stiffness = stiffness(shape_functions_derivatives)
reduced = np.delete(np.delete(stiffness, np.where(constraints), axis=0), np.where(constraints), axis=1)
loads_reduced = np.delete(loads, np.where(constraints))
displacements = np.linalg.solve(reduced, loads_reduced)
total_displacements = np.zeros(24)
total_displacements[np.where(constraints)] = 0
total_displacements[np.where(constraints == False)] = displacements

for idx, nid in enumerate(elem_ids):
    x, y, z = nodes[nid]
    # print(f'Node {nid}: {x, y, z}')
    print(f'u = {total_displacements[3 * idx]:.5f}, v = {total_displacements[3 * idx + 1]:.5f}, w = {total_displacements[3 * idx + 2]:.5f}')
