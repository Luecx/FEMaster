import sympy as sp
import numpy as np

# Define symbols for natural coordinates
g, h, r = sp.symbols('g h r')

# Create the nodes with shape functions and their locations
def node_definitions():
    nodes = [
        {"shape_function": -1/8 * (1 - g) * (1 - h) * (1 - r) * (2 + g + h + r), "location": (-1, -1, -1)},
        {"shape_function": -1/8 * (1 + g) * (1 - h) * (1 - r) * (2 - g + h + r), "location": (1, -1, -1)},
        {"shape_function": -1/8 * (1 + g) * (1 + h) * (1 - r) * (2 - g - h + r), "location": (1, 1, -1)},
        {"shape_function": -1/8 * (1 - g) * (1 + h) * (1 - r) * (2 + g - h + r), "location": (-1, 1, -1)},
        {"shape_function": -1/8 * (1 - g) * (1 - h) * (1 + r) * (2 + g + h - r), "location": (-1, -1, 1)},
        {"shape_function": -1/8 * (1 + g) * (1 - h) * (1 + r) * (2 - g + h - r), "location": (1, -1, 1)},
        {"shape_function": -1/8 * (1 + g) * (1 + h) * (1 + r) * (2 - g - h - r), "location": (1, 1, 1)},
        {"shape_function": -1/8 * (1 - g) * (1 + h) * (1 + r) * (2 + g - h - r), "location": (-1, 1, 1)},

        {"shape_function": 1/4 * (1 - g) * (1 + g) * (1 - h) * (1 - r), "location": (0, -1, -1)},
        {"shape_function": 1/4 * (1 - h) * (1 + h) * (1 + g) * (1 - r), "location": (1, 0, -1)},
        {"shape_function": 1/4 * (1 - g) * (1 + g) * (1 + h) * (1 - r), "location": (0, 1, -1)},
        {"shape_function": 1/4 * (1 - h) * (1 + h) * (1 - g) * (1 - r), "location": (-1, 0, -1)},
        {"shape_function": 1/4 * (1 - g) * (1 + g) * (1 - h) * (1 + r), "location": (0, -1, 1)},
        {"shape_function": 1/4 * (1 - h) * (1 + h) * (1 + g) * (1 + r), "location": (1, 0, 1)},
        {"shape_function": 1/4 * (1 - g) * (1 + g) * (1 + h) * (1 + r), "location": (0, 1, 1)},
        {"shape_function": 1/4 * (1 - h) * (1 + h) * (1 - g) * (1 + r), "location": (-1, 0, 1)},

        {"shape_function": 1/4 * (1 - r) * (1 + r) * (1 - g) * (1 - h), "location": (-1, -1, 0)},
        {"shape_function": 1/4 * (1 - r) * (1 + r) * (1 + g) * (1 - h), "location": (1, -1, 0)},
        {"shape_function": 1/4 * (1 - r) * (1 + r) * (1 + g) * (1 + h), "location": (1, 1, 0)},
        {"shape_function": 1/4 * (1 - r) * (1 + r) * (1 - g) * (1 + h), "location": (-1, 1, 0)}
    ]

    return nodes

def integration_scheme():
    # Gauss points for a 2x2x2 scheme (in natural coordinates, -1/sqrt(3), 1/sqrt(3))
    sqrt3 = np.sqrt(3) / 3
    gauss_points = [
        {"location": (-sqrt3, -sqrt3, -sqrt3), "weight": 1},
        {"location": ( sqrt3, -sqrt3, -sqrt3), "weight": 1},
        {"location": ( sqrt3,  sqrt3, -sqrt3), "weight": 1},
        {"location": (-sqrt3,  sqrt3, -sqrt3), "weight": 1},
        {"location": (-sqrt3, -sqrt3,  sqrt3), "weight": 1},
        {"location": ( sqrt3, -sqrt3,  sqrt3), "weight": 1},
        {"location": ( sqrt3,  sqrt3,  sqrt3), "weight": 1},
        {"location": (-sqrt3,  sqrt3,  sqrt3), "weight": 1}
    ]

    return gauss_points


def material_definitions(E, nu):
    # Calculate Lame parameters
    factor = E / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))

    # Create the 6x6 matrix using the Lame parameters
    C = np.array([
        [lambda_ + 2*mu, lambda_, lambda_, 0, 0, 0],
        [lambda_, lambda_ + 2*mu, lambda_, 0, 0, 0],
        [lambda_, lambda_, lambda_ + 2*mu, 0, 0, 0],
        [0, 0, 0, mu, 0, 0],
        [0, 0, 0, 0, mu, 0],
        [0, 0, 0, 0, 0, mu]
    ])

    # Return the matrix as a sparse matrix
    return sp.Matrix(C)


# Create a vector of shape functions
def shape_functions():
    nodes = node_definitions()
    # Extract shape functions from nodes and construct a column vector
    shape_function_vector = sp.Matrix([node['shape_function'] for node in nodes])
    return shape_function_vector

# Create the derivatives matrix (Nx3 matrix with derivatives w.r.t g, h, and r)
def local_derivative():
    funcs = shape_functions()

    # get all variables and determine if 2d or 3d
    variables = funcs.free_symbols

    # Initialize the matrix to store the derivatives (N x 3 matrix)
    derivatives_matrix = sp.zeros(len(funcs), len(variables))

    for i, func in enumerate(funcs):
        for j, variable in enumerate(variables):
            derivatives_matrix[i, j] = sp.diff(func, variable)

    return derivatives_matrix

# Create the strain-displacement matrix
def strain_displacement(derivatives_matrix):
    N = derivatives_matrix.rows
    C = derivatives_matrix.cols
    if C == 3:
        strain_displacement_matrix = sp.zeros(6, 3 * N)
        for i in range(N):
            strain_displacement_matrix[0, 3*i]   = derivatives_matrix[i, 0]  # dN/dg
            strain_displacement_matrix[1, 3*i+1] = derivatives_matrix[i, 1]  # dN/dh
            strain_displacement_matrix[2, 3*i+2] = derivatives_matrix[i, 2]  # dN/dr
            strain_displacement_matrix[3, 3*i]   = derivatives_matrix[i, 1]  # dN/dh
            strain_displacement_matrix[3, 3*i+1] = derivatives_matrix[i, 0]  # dN/dg
            strain_displacement_matrix[4, 3*i+1] = derivatives_matrix[i, 2]  # dN/dr
            strain_displacement_matrix[4, 3*i+2] = derivatives_matrix[i, 1]  # dN/dh
            strain_displacement_matrix[5, 3*i]   = derivatives_matrix[i, 2]  # dN/dr
            strain_displacement_matrix[5, 3*i+2] = derivatives_matrix[i, 0]  # dN/dg

    elif C == 2:
        strain_displacement_matrix = sp.zeros(3, 2 * N)
        for i in range(N):
            strain_displacement_matrix[0, 2*i]   = derivatives_matrix[i, 0]  # dN/dg
            strain_displacement_matrix[1, 2*i+1] = derivatives_matrix[i, 1]  # dN/dh
            strain_displacement_matrix[2, 2*i]   = derivatives_matrix[i, 1]  # dN/dh
            strain_displacement_matrix[2, 2*i+1] = derivatives_matrix[i, 0]  # dN/dg
    else:
        raise ValueError("Dimensionality must be 2 or 3")

    return strain_displacement_matrix


def compute_stiffness_matrix(E, nu):
    # Get the material matrix (constant)
    C = material_definitions(E, nu)

    # Initialize stiffness matrix K
    num_nodes = len(node_definitions())
    stiffness_matrix = sp.zeros(3 * num_nodes, 3 * num_nodes)

    # Get shape functions derivatives matrix
    derivatives_matrix = local_derivative()

    # Get the integration points
    gauss_points = integration_scheme()

    # Loop over all integration points
    for point in gauss_points:
        loc = point['location']
        weight = point['weight']

        # Substitute the gauss point location into the derivatives matrix
        B = strain_displacement(derivatives_matrix.subs({g: loc[0], h: loc[1], r: loc[2]}))

        # Compute the local stiffness contribution: K_local = B^T * C * B * det(Jacobian) * weight
        # Assuming det(Jacobian) = 1 for simplicity (unit cube element)
        K_local = B.T * C * B * weight

        # Sum this contribution to the global stiffness matrix
        stiffness_matrix += K_local

    return stiffness_matrix

# Testing the creation of the shape function vector and derivatives matrix
shape_function_vector = shape_functions()
derivatives_matrix = local_derivative()
strain_disp = strain_displacement(derivatives_matrix)


# # Display the results
# print("Shape Function Vector:")
# sp.pprint(shape_function_vector.subs({g: 0, h: 0, r: 0}))
# print("\nDerivatives Matrix (Nx3):")
# sp.pprint(derivatives_matrix.subs({g: 0, h: 0, r: 0}))
# print("\nStrain-Displacement Matrix:")
# sp.pprint(strain_disp.subs({g: 0, h: 0, r: 0}))

# material matrix
print(sp.pprint(material_definitions(1, 0)))

#
stiffness_matrix = compute_stiffness_matrix(1, 0)

# Print or pretty print the final stiffness matrix (symbolic)
print(sp.pprint(stiffness_matrix[0:6, 0:6] * 5))