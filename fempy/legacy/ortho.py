import numpy as np

def transformation_matrix(R):

    # Extract elements of R for readability
    R11, R12, R13 = R[0, 0], R[0, 1], R[0, 2]
    R21, R22, R23 = R[1, 0], R[1, 1], R[1, 2]
    R31, R32, R33 = R[2, 0], R[2, 1], R[2, 2]

    T_sigma = np.array([
        [R11**2, R21**2, R31**2, R11*R21, R21*R31, R31*R11],
        [R12**2, R22**2, R32**2, R12*R22, R22*R32, R32*R12],
        [R13**2, R23**2, R33**2, R13*R23, R23*R33, R33*R13],
        [R11*R12, R21*R22, R31*R32, R11*R22 + R12*R21, R21*R32 + R31*R22, R31*R12 + R32*R11],
        [R12*R13, R22*R23, R33*R32, R23*R12 + R13*R22, R22*R33 + R32*R23, R32*R13 + R12*R33],
        [R11*R13, R23*R21, R33*R31, R13*R21 + R11*R23, R23*R31 + R21*R33, R33*R11 + R31*R13]
    ])

    # Construct T_epsilon based on the specific pattern for epsilon transformation
    T_epsilon = T_sigma.copy()
    T_sigma  [:3, 3:] *= 2
    T_epsilon[3:, :3] *= 2

    return T_sigma, T_epsilon

def orthotropic_matrix(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23):
    # Compute auxiliary Poisson ratios for symmetry
    nu21 = nu12 * E2 / E1
    nu31 = nu13 * E3 / E1
    nu32 = nu23 * E3 / E2

    # Construct the compliance matrix S (inverse of the stiffness matrix C)
    S = np.zeros((6, 6))
    S[0, 0] = 1 / E1
    S[1, 1] = 1 / E2
    S[2, 2] = 1 / E3
    S[0, 1] = -nu12 / E1
    S[1, 0] = -nu21 / E2
    S[0, 2] = -nu13 / E1
    S[2, 0] = -nu31 / E3
    S[1, 2] = -nu23 / E2
    S[2, 1] = -nu32 / E3
    S[3, 3] = 1 / G23
    S[4, 4] = 1 / G13
    S[5, 5] = 1 / G12

    # Invert the compliance matrix to get the stiffness matrix
    C = np.linalg.inv(S)

    return C

def isotropic_matrix(E, nu):

    C = np.matrix([
        1-nu, nu, nu, 0, 0, 0,
        nu, 1-nu, nu, 0, 0, 0,
        nu, nu, 1-nu, 0, 0, 0,
        0, 0, 0, (1-2*nu)/2, 0, 0,
        0, 0, 0, 0, (1-2*nu)/2, 0,
        0, 0, 0, 0, 0, (1-2*nu)/2]
    ) * E / (1+nu) / (1-2*nu)
    return C.reshape((6, 6))


def rotation_matrix(t1, t2, t3):
    r1 = np.array([[np.cos(t1), np.sin(t1), 0]
                  , [-np.sin(t1), np.cos(t1), 0]
                  , [0, 0, 1]])
    r2 = np.array([[1, 0, 0]
                    , [0, np.cos(t2), np.sin(t2)]
                    , [0, -np.sin(t2), np.cos(t2)]])
    r3 = np.array([[np.cos(t3), np.sin(t3), 0]
                    , [-np.sin(t3), np.cos(t3), 0]
                    , [0, 0, 1]])
    return np.dot(np.dot(r1, r2), r3)

T_sigma, T_epsilon = transformation_matrix(rotation_matrix(35, 0, 0))
C    = orthotropic_matrix(210, 210, 30, 0.3, 0.3, 0.3, 210 / (2 * (1+0.3)), 210 / (2 * (1+0.3)), 210 / (2 * (1+0.3)))
# C2  = isotropic_matrix(210, 0.3)

np.set_printoptions(precision=2, suppress=True)

# print(M)
# print(N)
#
# print(M @ C @ M.T)
# print(N @ C @ N.T)

# print inverse(M) * C * N
print(T_epsilon.T @ C @ T_epsilon)

