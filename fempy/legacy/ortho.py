import numpy as np

def orthotropic_stiffness(E_x, E_y, nu_xy, G_xy):
    """
    Defines a 2x2x2x2 stiffness tensor for a 2D orthotropic material.

    Parameters:
    - E_x: Young's modulus in x direction
    - E_y: Young's modulus in y direction
    - nu_xy: Poisson's ratio in the xy plane
    - G_xy: Shear modulus in the xy plane

    Returns:
    - C: The 2x2x2x2 stiffness tensor for orthotropic material
    """
    nu_yx = nu_xy * E_y / E_x
    C = np.zeros((2, 2, 2, 2))
    C[0, 0, 0, 0] = E_x / (1 - nu_xy * nu_yx)
    C[1, 1, 1, 1] = E_y / (1 - nu_xy * nu_yx)
    C[0, 0, 1, 1] = C[1, 1, 0, 0] = nu_xy * E_y / (1 - nu_xy * nu_yx)
    C[0, 1, 0, 1] = C[1, 0, 0, 1] = C[0, 1, 1, 0] = C[1, 0, 1, 0] = G_xy
    return C

def stiffness_matrix_3x3(E_x, E_y, nu_xy, G_xy):
    """
    Converts the 2x2x2x2 stiffness tensor to a 3x3 stiffness matrix for 2D orthotropic materials.

    Parameters:
    - E_x: Young's modulus in x direction
    - E_y: Young's modulus in y direction
    - nu_xy: Poisson's ratio in the xy plane
    - G_xy: Shear modulus in the xy plane

    Returns:
    - C_mat: The 3x3 stiffness matrix for orthotropic material
    """
    nu_yx = nu_xy * E_y / E_x

    # Define components in the 3x3 stiffness matrix form
    C11 = E_x / (1 - nu_xy * nu_yx)
    C22 = E_y / (1 - nu_xy * nu_yx)
    C12 = nu_xy * E_y / (1 - nu_xy * nu_yx)
    C33 = G_xy

    # Construct the 3x3 matrix
    C_mat = np.array([
        [C11, C12, 0],
        [C12, C22, 0],
        [0, 0, C33]
    ])

    return C_mat

def strain_vector_and_tensor(e_xx, e_yy, gamma_xy):
    """
    Returns the strain vector and strain tensor for given strain parameters in 2D.

    Parameters:
    - e_xx: Normal strain in the x-direction
    - e_yy: Normal strain in the y-direction
    - gamma_xy: Shear strain in the xy-plane

    Returns:
    - strain_vector: The strain vector [e_xx, e_yy, gamma_xy]
    - strain_tensor: The strain tensor in 2x2 matrix form
    """
    # Define the strain vector
    strain_vector = np.array([e_xx, e_yy, 2 * gamma_xy])

    # Define the strain tensor (note gamma_xy / 2 for the off-diagonal terms)
    strain_tensor = np.array([
        [e_xx, gamma_xy ],
        [gamma_xy, e_yy]
    ])

    return strain_vector, strain_tensor

# Example usage
E_x = 150.0   # Young's modulus in x direction (MPa)
E_y = 120.0   # Young's modulus in y direction (MPa)
nu_xy = 0.25  # Poisson's ratio in the xy plane
G_xy = 50.0   # Shear modulus in xy plane (MPa)

# Generate stiffness tensor and matrix
C_tensor = orthotropic_stiffness(E_x, E_y, nu_xy, G_xy)
C_mat = stiffness_matrix_3x3(E_x, E_y, nu_xy, G_xy)

# print("2x2x2x2 stiffness tensor for orthotropic material:\n", C_tensor)
# print("\n3x3 stiffness matrix for orthotropic material:\n", C_mat)

# Define strain parameters
e_xx = 0.01   # Strain in x direction
e_yy = 0.005  # Strain in y direction
gamma_xy = 0.002  # Shear strain in xy plane

# Generate strain vector and tensor
strain_vec, strain_tensor = strain_vector_and_tensor(e_xx, e_yy, gamma_xy)
# print("\nStrain vector:\n", strain_vec)
# print("\nStrain tensor:\n", strain_tensor)

print(C_mat)
# print(C_tensor)

# stress tensor and vector
stress_vec = np.dot(C_mat, strain_vec)
stress_tensor = np.tensordot(C_tensor, strain_tensor, axes=([2, 3], [0, 1]))

print("\nStress vector:\n", stress_vec)
print("\nStress tensor:\n", stress_tensor)

