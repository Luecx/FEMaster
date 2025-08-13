import numpy as np

class S4:
    def __init__(self, node_coords):
        """
        Initializes the shell element with 4 nodes.

        Parameters:
        node_coords (list of np.array): List of 4 node coordinates, each a 3D vector.
        """
        assert len(node_coords) == 4, "S4 element requires exactly 4 nodes."
        self.node_coords = np.array(node_coords)
        self.local_coords = self.local_coordinates()

    def shape_function(self, r, s):
        """
        Computes the shape function values at local coordinates (r, s).

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: Shape function values.
        """
        return np.array([
            (1 - r) * (1 - s) / 4,
            (1 + r) * (1 - s) / 4,
            (1 + r) * (1 + s) / 4,
            (1 - r) * (1 + s) / 4,
            ])

    def shape_function_derivative(self, r, s):
        """
        Computes the derivatives of shape functions with respect to (r, s).

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: Derivatives of shape functions with respect to r and s.
        """
        return np.array([
            [-(1 - s) / 4, -(1 - r) / 4],
            [ (1 - s) / 4, -(1 + r) / 4],
            [ (1 + s) / 4,  (1 + r) / 4],
            [-(1 + s) / 4,  (1 - r) / 4],
        ])

    def jacobian(self, r, s):
        """
        Computes the Jacobian matrix at local coordinates (r, s).

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: 2x2 Jacobian matrix.
        """
        dN = self.shape_function_derivative(r, s)
        J = dN.T @ self.local_coords[:, :2]
        return J

    def strain_displacement_bending(self, r, s):
        """
        Computes the strain-displacement matrix for bending.

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: Strain-displacement matrix for bending.
        """
        dN = self.shape_function_derivative(r, s)
        J = self.jacobian(r, s)
        J_inv = np.linalg.inv(J)
        dH = (dN @ J_inv.T).T

        print("evaluating at r,s = ", r, s)
        print("dN = ", dN)
        print("J = ", J)
        print("J_inv = ", J_inv)

        B_bend = np.zeros((3, 12))
        for i in range(4):
            B_bend[0, 3 * i + 2] = -dH[0, i]
            B_bend[1, 3 * i + 1] =  dH[1, i]
            B_bend[2, 3 * i + 1] =  dH[0, i]
            B_bend[2, 3 * i + 2] = -dH[1, i]

        return B_bend

    def strain_displacement_shear(self, r, s):
        """
        Computes the strain-displacement matrix for shear.

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: Strain-displacement matrix for shear.
        """
        dN = self.shape_function_derivative(r, s)
        J = self.jacobian(r, s)
        J_inv = np.linalg.inv(J)
        dH = (dN @ J_inv.T).T
        H = self.shape_function(r, s)

        B_shear = np.zeros((2, 12))
        for i in range(4):
            B_shear[0, 3 * i] = dH[0, i]
            B_shear[0, 3 * i + 2] = H[i]
            B_shear[1, 3 * i] = dH[1, i]
            B_shear[1, 3 * i + 1] = -H[i]

        return B_shear

    def strain_displacement_membrane(self, r, s):
        """
        Computes the strain-displacement matrix for membrane.

        Parameters:
        r, s (float): Local coordinates in the reference domain.

        Returns:
        np.array: Strain-displacement matrix for membrane.
        """
        dN = self.shape_function_derivative(r, s)
        J = self.jacobian(r, s)
        J_inv = np.linalg.inv(J)
        dH = (dN @ J_inv.T).T

        B_membrane = np.zeros((3, 8))
        for i in range(4):
            B_membrane[0, 2 * i] = dH[0, i]
            B_membrane[1, 2 * i + 1] = dH[1, i]
            B_membrane[2, 2 * i] = dH[1, i]
            B_membrane[2, 2 * i + 1] = dH[0, i]

        return B_membrane

    def stiffness_bending(self, elasticity_matrix):
        """
        Computes the bending stiffness matrix.

        Parameters:
        elasticity_matrix (np.array): Elasticity matrix for bending.
        thickness (float): Thickness of the shell.

        Returns:
        np.array: Stiffness matrix for bending.
        """
        gauss_points = [
            (-1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3),  1 / np.sqrt(3)),
            (-1 / np.sqrt(3),  1 / np.sqrt(3)),
        ]
        weights = [1, 1, 1, 1]

        K_bend = np.zeros((12, 12))
        for (r, s), w in zip(gauss_points, weights):
            J = self.jacobian(r, s)
            det_J = np.linalg.det(J)
            B_bend = self.strain_displacement_bending(r, s)

            K_bend += w * B_bend.T @ elasticity_matrix @ B_bend * det_J
        # bending correlates


        indices = [2, 3, 4,
               2+6, 3+6, 4+6,
               2+12, 3+12, 4+12,
               2+18, 3+18, 4+18]

        K_bend_loc = np.zeros((24, 24))
        K_bend_loc[np.ix_(indices, indices)] = K_bend

        return K_bend_loc

    def stiffness_shear(self, elasticity_matrix):
        """
        Computes the shear stiffness matrix.

        Parameters:
        elasticity_matrix (np.array): Elasticity matrix for shear.
        thickness (float): Thickness of the shell.

        Returns:
        np.array: Stiffness matrix for shear.
        """
        gauss_points = [
            (-1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3),  1 / np.sqrt(3)),
            (-1 / np.sqrt(3),  1 / np.sqrt(3)),
        ]
        weights = [1, 1, 1, 1]

        K_shear = np.zeros((12, 12))
        for (r, s), w in zip(gauss_points, weights):
            J = self.jacobian(r, s)
            det_J = np.linalg.det(J)
            B_shear = self.strain_displacement_shear(r, s)

            K_shear += w * B_shear.T @ elasticity_matrix @ B_shear * det_J

        indices = [2, 3, 4,
                   2+6, 3+6, 4+6,
                   2+12, 3+12, 4+12,
                   2+18, 3+18, 4+18]

        K_bend_loc = np.zeros((24, 24))
        K_bend_loc[np.ix_(indices, indices)] = K_shear

        return K_bend_loc

    def stiffness_membrane(self, elasticity_matrix):
        """
        Computes the membrane stiffness matrix.

        Parameters:
        elasticity_matrix (np.array): Elasticity matrix for membrane.
        thickness (float): Thickness of the shell.

        Returns:
        np.array: Stiffness matrix for membrane.
        """
        gauss_points = [
            (-1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3), -1 / np.sqrt(3)),
            ( 1 / np.sqrt(3),  1 / np.sqrt(3)),
            (-1 / np.sqrt(3),  1 / np.sqrt(3)),
        ]
        weights = [1, 1, 1, 1]


        K_membrane = np.zeros((8, 8))
        for (r, s), w in zip(gauss_points, weights):
            J = self.jacobian(r, s)
            det_J = np.linalg.det(J)
            B_membrane = self.strain_displacement_membrane(r, s)

            K_membrane += w * B_membrane.T @ elasticity_matrix @ B_membrane * det_J

        indices = [0, 1,
                   6, 7,
                   12, 13,
                   18, 19]
        K_membrane_loc = np.zeros((24, 24))
        K_membrane_loc[np.ix_(indices, indices)] = K_membrane
        return K_membrane_loc

    def area(self):
        gp = 1/np.sqrt(3)
        J1 = np.linalg.det(self.jacobian(-gp, -gp))
        J2 = np.linalg.det(self.jacobian( gp, -gp))
        J3 = np.linalg.det(self.jacobian( gp,  gp))
        J4 = np.linalg.det(self.jacobian(-gp,  gp))
        return (J1 + J2 + J3 + J4)

    def transformation(self):
        """
        Computes the transformation matrix for the element.

        Returns:
        np.array: Transformation matrix.
        """
        n1, n2, n3, n4 = self.node_coords[:4]
        n12 = n2 - n1
        n13 = n3 - n1

        x_axis = n12 / np.linalg.norm(n12)
        z_axis = np.cross(n12, n13) / np.linalg.norm(np.cross(n12, n13))
        y_axis = np.cross(z_axis, x_axis)

        axes = np.vstack([x_axis, y_axis, z_axis])

        T = np.zeros((24, 24))
        for i in range(8):
            T[3 * i:3 * i + 3, 3 * i:3 * i + 3] = axes

        return T

    def local_coordinates(self):
        n1, n2, n3, n4 = self.node_coords[:4]
        n11 = n1 - n1
        n12 = n2 - n1
        n13 = n3 - n1
        n14 = n4 - n1

        x_axis = n12 / np.linalg.norm(n12)
        z_axis = np.cross(n12, n13) / np.linalg.norm(np.cross(n12, n13))
        y_axis = np.cross(z_axis, x_axis)

        n1_loc = np.dot(n11, x_axis), np.dot(n11, y_axis), np.dot(n11, z_axis)
        n2_loc = np.dot(n12, x_axis), np.dot(n12, y_axis), np.dot(n12, z_axis)
        n3_loc = np.dot(n13, x_axis), np.dot(n13, y_axis), np.dot(n13, z_axis)
        n4_loc = np.dot(n14, x_axis), np.dot(n14, y_axis), np.dot(n14, z_axis)

        # n1_loc = np.array([0.0, 0.0,0.0])
        # n2_loc = np.array([1.5297058540778354, 0.0,0.0])
        # n3_loc = np.array([0.5491251783869153, 1.8434916702989301,0.0])
        # n4_loc = np.array([-0.6667948594698258, 0.7844645405527361,0.0])

        return np.array([n1_loc, n2_loc, n3_loc, n4_loc])


    def total_stiffness(self, elasticity_matrices):
        """
        Computes the total stiffness matrix including bending, shear, and membrane.

        Parameters:
        elasticity_matrices (dict): Dictionary with keys 'bending', 'shear', and 'membrane'.
        thickness (float): Thickness of the shell.

        Returns:
        np.array: Total stiffness matrix.
        """
        K_bend = self.stiffness_bending(elasticity_matrices['bending'])
        K_shear = self.stiffness_shear(elasticity_matrices['shear'])
        K_membrane = self.stiffness_membrane(elasticity_matrices['membrane'])

        print(K_bend)
        exit(0)

        # K_local = K_bend + K_shear + K_membrane
        #
        # T = self.transformation()
        # K_total = T.T @ K_local @ T

        return K_total

# # Example usage
node_positions = [
    np.array([0.1,-0.5, 0.0]),
    np.array([1.5, 0.3, 0.0]),
    np.array([1.0, 1.2, 0.0]),
    np.array([-0.5, 1.5, 0.0]),
]

# # Example usage
# node_positions = [
#     np.array([0.0, 0.0,0.0]),
#     np.array([1.5297058540778354, 0.0,0.0]),
#     np.array([0.5491251783869153, 1.8434916702989301,0.0]),
#     np.array([-0.6667948594698258, 0.7844645405527361,0.0]),
# ]



E = 1  # Young's modulus (Pa)
nu = 0  # Poisson's ratio
h = 3   # Thickness (m)

elasticity_matrices = {
    'bending': E * h**3 / (12 * (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ]),
    'shear': E * h * (5 / 6) / (2 * (1 + nu)) * np.array([
        [1, 0],
        [0, 1]
    ]),
    'membrane': h*(1 / (1 - nu**2)) * np.array([
        [E, nu * E, 0],
        [nu * E, E, 0],
        [0, 0, E / (2 * (1 + nu))]
    ])
}


print(elasticity_matrices)

np.set_printoptions(precision=2, suppress=True, linewidth=10000)

shell = S4(node_positions)


K_total = shell.total_stiffness(elasticity_matrices)
print(K_total)
