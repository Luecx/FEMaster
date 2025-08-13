# import numpy as np
#
# E = 1
# A = 1
# L = 1
# I_y = 1
# I_z = 1
# J = 1
# G = 1
#
#
# K11 = np.matrix([[ E*A/L,        0,             0,           0,              0,             0         ],
#                  [ 0,      12*E*I_z/L**3,       0,           0,              0,    -6*E*I_z/L**2    ],
#                  [ 0,            0,      12*E*I_y/L**3,      0,     6*E*I_y/L**2,            0         ],
#                  [ 0,            0,             0,      G*J/L,              0,             0         ],
#                  [ 0,            0,      6*E*I_y/L**2,       0,      4*E*I_y/L,              0         ],
#                  [ 0,     -6*E*I_z/L**2,        0,           0,              0,       4*E*I_z/L     ]])
#
# K12 = np.matrix([[-E*A/L,        0,             0,           0,              0,             0         ],
#                  [ 0,     -12*E*I_z/L**3,       0,           0,              0,    -6*E*I_z/L**2     ],
#                  [ 0,            0,     -12*E*I_y/L**3,      0,     6*E*I_y/L**2,            0         ],
#                  [ 0,            0,             0,     -G*J/L,              0,             0         ],
#                  [ 0,            0,      6*E*I_y/L**2,       0,      2*E*I_y/L,              0         ],
#                  [ 0,     -6*E*I_z/L**2,        0,           0,              0,       2*E*I_z/L      ]])
#
#
# K = np.bmat([[K11, K12],
#              [K12.T, K11]])
#
# # constraint left node --> remove first 6 rows and cols
# K = K11
#
#
# for idx in range(6):
#     F = np.zeros((6, 1))
#     F[idx] = 1
#     u = np.linalg.solve(K, F)
#     print(u.T)

import numpy as np

# Define constants
E = 1   # Young's modulus
A = 1   # Cross-sectional area
L = 1   # Length of the beam
I_y = 1 # Moment of inertia around the y-axis
I_z = 1 # Moment of inertia around the z-axis
J = 1   # Polar moment of inertia for torsion
G = 1   # Shear modulus

# Number of integration points (for Gaussian quadrature)
n_gauss = 2
gauss_points = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])  # Gauss points
gauss_weights = np.array([1, 1])  # Weights for Gaussian quadrature


def shape_functions(xi):
    """
    Shape functions for the beam element in local coordinates (xi).
    xi is the local coordinate ranging from -1 to 1.
    """
    # Translational shape functions
    N1 = (1 - xi) / 2  # Shape function for node 1
    N2 = (1 + xi) / 2  # Shape function for node 2

    return np.array([N1, N2])


def bending_shape_functions(xi):
    """
    Bending shape functions for the beam element.
    xi is the local coordinate ranging from -1 to 1.
    These are Hermite shape functions for bending.
    """
    L_half = L / 2
    N1 = 1 - 3 * xi ** 2 + 2 * xi ** 3
    N2 = L_half * (xi - 2 * xi ** 2 + xi ** 3)
    N3 = 3 * xi ** 2 - 2 * xi ** 3
    N4 = L_half * (-xi ** 2 + xi ** 3)

    return np.array([N1, N2, N3, N4])


def strain_displacement_matrix(xi):
    """
    Strain-displacement matrix B for axial deformation (N').
    xi is the local coordinate ranging from -1 to 1.
    """
    # Derivatives of the shape functions
    dN1_dxi = -0.5
    dN2_dxi = 0.5

    B = np.array([dN1_dxi, dN2_dxi]) * 2 / L  # Convert to physical coordinates

    return B


def bending_B_matrix(xi):
    """
    B matrix for bending strain (second derivative of shape functions).
    """
    L_half = L / 2
    d2N1_dxi2 = -6 / L**2 + 12 * xi / L**2
    d2N2_dxi2 = 1 - 2 * xi
    d2N3_dxi2 = 6 / L**2 - 12 * xi / L**2
    d2N4_dxi2 = -xi + xi ** 2

    return np.array([d2N1_dxi2, d2N2_dxi2, d2N3_dxi2, d2N4_dxi2])


def stiffness_matrix():
    """
    Compute the stiffness matrix for the beam element using Gaussian quadrature.
    """
    K = np.zeros((6, 6))  # Stiffness matrix (6x6 for a beam element)

    # Material constitutive matrices for axial and bending
    D_axial = E * A
    D_bending_y = E * I_y
    D_bending_z = E * I_z
    D_torsion = G * J

    # Integration over the element using Gaussian quadrature
    for i in range(n_gauss):
        xi = gauss_points[i]
        w = gauss_weights[i]

        # Compute B matrices
        B_axial = strain_displacement_matrix(xi)
        B_bending_y = bending_B_matrix(xi)  # B matrix for bending about y
        B_bending_z = bending_B_matrix(xi)  # B matrix for bending about z

        # Axial stiffness contribution
        K_axial = w * B_axial.T * D_axial @ B_axial * L / 2
        K[:2, :2] += K_axial

        # Bending stiffness contribution (around y-axis and z-axis)
        K_bending_y = w * B_bending_y.T * D_bending_y @ B_bending_y * L / 2
        K_bending_z = w * B_bending_z.T * D_bending_z @ B_bending_z * L / 2
        K[2:6, 2:6] += K_bending_y + K_bending_z

    return K


# Compute the stiffness matrix using interpolation functions and numerical integration
K = stiffness_matrix()

print(K)

# # Example to solve for displacements with constrained left node (remove rows/cols)
# K_reduced = K[2:, 2:]  # Remove first 2 rows/cols to simulate the constraint
#
# # Solve for displacement for each force case
# for idx in range(4):
#     F = np.zeros((4, 1))
#     F[idx] = 1
#     u = np.linalg.solve(K_reduced, F)
#     print(u.T)
