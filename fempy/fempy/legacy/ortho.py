import numpy as np

def rotation_matrix(phi1, phi2, Phi):
    R = np.zeros((3, 3))
    R[0, 0] = np.cos(phi1) * np.cos(phi2) - np.cos(Phi) * np.sin(phi1) * np.sin(phi2)
    R[0, 1] = np.cos(phi1) * np.sin(phi2) + np.cos(Phi) * np.sin(phi1) * np.cos(phi2)
    R[0, 2] = np.sin(Phi) * np.sin(phi1)
    R[1, 0] = -np.sin(phi1) * np.cos(phi2) - np.cos(Phi) * np.cos(phi1) * np.sin(phi2)
    R[1, 1] = -np.sin(phi1) * np.sin(phi2) + np.cos(Phi) * np.cos(phi1) * np.cos(phi2)
    R[1, 2] = np.sin(Phi) * np.cos(phi1)
    R[2, 0] = np.sin(Phi) * np.sin(phi2)
    R[2, 1] = -np.sin(Phi) * np.cos(phi2)
    R[2, 2] = np.cos(Phi)
    return R

def stiffness_matrix():
    E = 230
    nu = 0.3



    C = np.array([
        [1-nu, nu, nu, 0, 0, 0],
        [nu, 1-nu, nu, 0, 0, 0],
        [nu, nu, 1-nu, 0, 0, 0],
        [0, 0, 0, (1-2*nu) / 2, 0, 0],
        [0, 0, 0, 0, (1-2*nu) / 2, 0],
        [0, 0, 0, 0, 0, (1-2*nu) / 2]
    ])   # Convert GPa to Pa
    C *= E / (1 + nu) / (1 - 2 * nu)
    return C

def transformation_matrix(R):
    # Define T_sigma and T_epsilon here based on the components of R
    T_sigma = np.zeros((6, 6))
    T_epsilon = np.zeros((6, 6))

    # Fill in T_sigma and T_epsilon according to the formulas in the image
    # Note: This part requires writing each component based on R
    # Assigning values to T_sigma
    T_sigma[0, 0] = R[0, 0] ** 2
    T_sigma[0, 1] = R[1, 0] ** 2
    T_sigma[0, 2] = R[2, 0] ** 2
    T_sigma[0, 3] = R[0, 0] * R[1, 0]
    T_sigma[0, 4] = R[1, 0] * R[2, 0]
    T_sigma[0, 5] = R[2, 0] * R[0, 0]

    T_sigma[1, 0] = R[0, 1] ** 2
    T_sigma[1, 1] = R[1, 1] ** 2
    T_sigma[1, 2] = R[2, 1] ** 2
    T_sigma[1, 3] = R[0, 1] * R[1, 1]
    T_sigma[1, 4] = R[1, 1] * R[2, 1]
    T_sigma[1, 5] = R[2, 1] * R[0, 1]

    T_sigma[2, 0] = R[0, 2] ** 2
    T_sigma[2, 1] = R[1, 2] ** 2
    T_sigma[2, 2] = R[2, 2] ** 2
    T_sigma[2, 3] = R[0, 2] * R[1, 2]
    T_sigma[2, 4] = R[1, 2] * R[2, 2]
    T_sigma[2, 5] = R[2, 2] * R[0, 2]

    T_sigma[3, 0] = R[0, 0] * R[0, 1]
    T_sigma[3, 1] = R[1, 0] * R[1, 1]
    T_sigma[3, 2] = R[2, 0] * R[2, 1]
    T_sigma[3, 3] = R[0, 0] * R[1, 1] + R[1, 0] * R[0, 1]
    T_sigma[3, 4] = R[1, 1] * R[2, 0] + R[2, 1] * R[1, 0]
    T_sigma[3, 5] = R[2, 0] * R[0, 1] + R[0, 0] * R[2, 1]

    T_sigma[4, 0] = R[0, 1] * R[0, 2]
    T_sigma[4, 1] = R[1, 1] * R[1, 2]
    T_sigma[4, 2] = R[2, 1] * R[2, 2]
    T_sigma[4, 3] = R[0, 1] * R[1, 2] + R[1, 1] * R[0, 2]
    T_sigma[4, 4] = R[1, 2] * R[2, 1] + R[2, 2] * R[1, 1]
    T_sigma[4, 5] = R[2, 1] * R[0, 2] + R[0, 1] * R[2, 2]

    T_sigma[5, 0] = R[0, 2] * R[0, 0]
    T_sigma[5, 1] = R[1, 2] * R[1, 0]
    T_sigma[5, 2] = R[2, 2] * R[2, 0]
    T_sigma[5, 3] = R[0, 2] * R[1, 0] + R[1, 2] * R[0, 0]
    T_sigma[5, 4] = R[1, 0] * R[2, 2] + R[2, 0] * R[1, 2]
    T_sigma[5, 5] = R[2, 0] * R[0, 2] + R[0, 0] * R[2, 2]

    T_epsilon = T_sigma.copy()

    T_epsilon[3:, :3] *= 2
    T_sigma[:3, 3:] *= 2

    return T_sigma, T_epsilon

def rotated_stiffness_matrix(C, T_sigma, T_epsilon):
    return np.dot(np.dot(np.linalg.inv(T_sigma), C), T_epsilon)

# Example usage:
phi1, phi2, Phi = np.radians([30, 156, 2])  # example angles in radians
R = rotation_matrix(phi1, phi2, Phi)

C = stiffness_matrix()

print(C)
print(R)

C = stiffness_matrix()
T_sigma, T_epsilon = transformation_matrix(R)
C_EBSP0 = rotated_stiffness_matrix(C, T_sigma, T_epsilon)


np.set_printoptions(precision=2,  suppress=True)


print("Rotated stiffness matrix C_EBSP0:\n", C_EBSP0)
