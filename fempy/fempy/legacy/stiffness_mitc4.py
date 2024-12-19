# copied from pynite to verify the correctness of the stiffness matrix
# adjusted to not use MITC but the standard B matrix

# The MITC4 element is no longer used by Pynite.
# It has been replaced with the DKMQ element instead. Still, it seems a waste to throw this legacy code out.
# There may be a use for an MITC4 element in the future - even if only for academic purposes.
# Many FEA programs allow the user to pick from multiple types of plate elements. Pynite may take that approach in the future.


# References used to derive this element:
# 1. "Finite Element Procedures, 2nd Edition", Klaus-Jurgen Bathe
# 2. "A First Course in the Finite Element Method, 4th Edition", Daryl L. Logan
# 3. "Finite Element Analysis Fundamentals", Richard H. Gallagher

from numpy import array, arccos, dot, cross, matmul, add, zeros
from numpy.linalg import inv, det, norm
from math import sin, cos

class Node():
    X = 0
    Y = 0
    Z = 0

    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z

class Material():
    E = 0
    nu = 0

    def __init__(self, E, nu):
        self.E = E
        self.nu = nu

class Model():
    nodes = {}
    materials = {}

    def __init__(self):

        a = 2

        self.nodes[1] = Node(0.6, 0.1, 0)
        self.nodes[2] = Node(2, -0.1, 0)
        self.nodes[3] = Node(2.3, 2, 0)
        self.nodes[4] = Node(-0.5, 2, 0)

        self.materials['Steel'] = Material(1, 0)

class MITC4():
    """
    An isoparametric general quadrilateral element, formulated by superimposing an isoparametric
    MITC4 bending element with an isoparametric plane stress element. Drilling stability is
    provided by adding a weak rotational spring stiffness at each node. Isotropic behavior is the
    default, but orthotropic in-plane behavior can be modeled by specifying stiffness modification
    factors for the element's local x and y axes.

    This element performs well for thick and thin plates, and for skewed plates. Element center
    stresses and corner FORCES converge rapidly; however, corner STRESSES are more representative
    of center stresses. Minor errors are introduced into the solution due to the drilling
    approximation. Orthotropic behavior is limited to acting along the plate's local axes.
    """

    #%%
    def __init__(self, name, i_node, j_node, m_node, n_node, t, material_name, model, kx_mod=1.0,
                 ky_mod=1.0):

        self.name = name
        self.ID = None
        self.type = 'MITC4'

        self.i_node = i_node
        self.j_node = j_node
        self.m_node = m_node
        self.n_node = n_node

        self.t = t
        self.kx_mod = kx_mod
        self.ky_mod = ky_mod

        self.pressures = []  # A list of surface pressures [pressure, case='Case 1']

        # Quads need a link to the model they belong to
        self.model = model

        # Get material properties for the plate from the model
        try:
            self.E = self.model.materials[material_name].E
            self.nu = self.model.materials[material_name].nu
        except:
            raise KeyError('Please define the material ' + str(material_name) + ' before assigning it to plates.')

    #%%
    def _local_coords(self):
        '''
        Calculates or recalculates and stores the local (x, y) coordinates for each node of the
        quadrilateral.
        '''

        # Get the global coordinates for each node
        X1, Y1, Z1 = self.m_node.X, self.m_node.Y, self.m_node.Z
        X2, Y2, Z2 = self.n_node.X, self.n_node.Y, self.n_node.Z
        X3, Y3, Z3 = self.i_node.X, self.i_node.Y, self.i_node.Z
        X4, Y4, Z4 = self.j_node.X, self.j_node.Y, self.j_node.Z


        # Following Reference 1, Figure 5.26, node 3 will be used as the
        # origin of the plate's local (x, y) coordinate system. Find the
        # vector from the origin to each node.
        vector_32 = array([X2 - X3, Y2 - Y3, Z2 - Z3]).T
        vector_31 = array([X1 - X3, Y1 - Y3, Z1 - Z3]).T
        vector_34 = array([X4 - X3, Y4 - Y3, Z4 - Z3]).T

        # Define the plate's local x, y, and z axes
        x_axis = vector_34
        z_axis = cross(x_axis, vector_32)
        y_axis = cross(z_axis, x_axis)

        # Convert the x and y axes into unit vectors
        x_axis = x_axis/norm(x_axis)
        y_axis = y_axis/norm(y_axis)

        # Calculate the local (x, y) coordinates for each node
        self.x1 = dot(vector_31, x_axis)
        self.x2 = dot(vector_32, x_axis)
        self.x3 = 0
        self.x4 = dot(vector_34, x_axis)
        self.y1 = dot(vector_31, y_axis)
        self.y2 = dot(vector_32, y_axis)
        self.y3 = 0
        self.y4 = dot(vector_34, y_axis)

    #%%
    def J(self, r, s):
        '''
        Returns the Jacobian matrix for the element
        '''

        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4

        # Return the Jacobian matrix
        return 1/4*array([[x1*(s + 1) - x2*(s + 1) + x3*(s - 1) - x4*(s - 1), y1*(s + 1) - y2*(s + 1) + y3*(s - 1) - y4*(s - 1)],
                          [x1*(r + 1) - x2*(r - 1) + x3*(r - 1) - x4*(r + 1), y1*(r + 1) - y2*(r - 1) + y3*(r - 1) - y4*(r + 1)]])

    #%%
    def B_kappa(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[1 + s, -1 - s, -1 + s,  1 - s],
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))

        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_kappa = array([[0,    0,     -dH[0, 0], 0,    0,     -dH[0, 1], 0,    0,     -dH[0, 2], 0,    0,     -dH[0, 3]],
                         [0, dH[1, 0],     0,     0, dH[1, 1],     0,     0, dH[1, 2],     0,     0, dH[1, 3],     0    ],
                         [0, dH[0, 0], -dH[1, 0], 0, dH[0, 1], -dH[1, 1], 0, dH[0, 2], -dH[1, 2], 0, dH[0, 3], -dH[1, 3]]])

        # Below is the matrix derived from the 1984 version of the MITC4 element. It appears to be
        # the same, but with a different sign convention for the section curvatures.
        # B_kappa = array([[0,     0,     dH[0, 0],  0,     0,     dH[0, 1],  0,     0,     dH[0, 2],  0,     0,     dH[0, 3]],
        #                  [0,  dH[1, 0],     0,     0,  dH[1, 1],     0,     0,  dH[1, 2],     0,     0,  dH[1, 3],     0   ],
        #                  [0, -dH[0, 0], dH[1, 0],  0, -dH[0, 1], dH[1, 1],  0, -dH[0, 2], dH[1, 2],  0, -dH[0, 3], dH[1, 3]]])

        return B_kappa

    #%%
    def B_gamma(self, r, s):
        '''
        Returns the [B] matrix for shear.

        This is provided for reference only and is not actually used by
        PyNite. This is the theoretical solution, but it is known to
        produce spurious shear forces. It is prone to a phenomenon called
        shear locking. Instead of this matrix, the MITC4 [B] matrix is used,
        which eliminates shear-locking and can be used for thick and thin
        plates.
        '''

        H = 1/4*array([(1 + r)*(1 + s), (1 - r)*(1 + s), (1 - r)*(1 - s), (1 + r)*(1 - s)])

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with respect to r and s
        # to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[1 + s, -1 - s, -1 + s,  1 - s],
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))

        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_gamma = array([[dH[0, 0],   0,   H[0], dH[0, 1],   0,   H[1], dH[0, 2],   0,   H[2], dH[0, 3],   0,   H[3]],
                         [dH[1, 0], -H[0],  0,   dH[1, 1], -H[1],  0,   dH[1, 2], -H[2],  0,   dH[1, 3], -H[3],  0  ]])


        return B_gamma

    def B_gamma_MITC4(self, r, s):
        '''
        Returns the [B] matrix for shear.

        MITC stands for mixed interpolation tensoral components. MITC elements
        are used in many programs and are known to perform well for thick and
        thin plates, and for distorted plate geometries.
        '''

        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4
        x_axis = array([1, 0, 0]).T

        # Reference 1, Equations 5.105
        Ax = x1 - x2 - x3 + x4
        Bx = x1 - x2 + x3 - x4
        Cx = x1 + x2 - x3 - x4
        Ay = y1 - y2 - y3 + y4
        By = y1 - y2 + y3 - y4
        Cy = y1 + y2 - y3 - y4

        # Find the angles between the axes of the natural coordinate system and
        # the local x-axis.
        r_axis = array([(x1 + x4)/2 - (x2 + x3)/2, (y1 + y4)/2 - (y2 + y3)/2, 0]).T
        s_axis = array([(x1 + x2)/2 - (x3 + x4)/2, (y1 + y2)/2 - (y3 + y4)/2, 0]).T

        r_axis = r_axis/norm(r_axis)
        s_axis = s_axis/norm(s_axis)

        alpha = arccos(dot(r_axis, x_axis))
        beta = arccos(dot(s_axis, x_axis))
        # alpha = atan(Ay/Ax)
        # beta = pi/2 - atan(Cx/Cy)

        # Reference 1, Equations 5.103 and 5.104 (p. 426)
        det_J = det(self.J(r, s))

        gr = ((Cx + r*Bx)**2 + (Cy + r*By)**2)**0.5/(8*det_J)
        gs = ((Ax + s*Bx)**2 + (Ay + s*By)**2)**0.5/(8*det_J)

        # d      =           [    w1           theta_x1             theta_y1             w2            theta_x2              theta_y2            w3             theta_x3             theta_y3         w4             theta_x4             theta_y4      ]
        gamma_rz = gr*array([[(1 + s)/2, -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 + s)/2,  -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 - s)/2, -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s), (1 - s)/2,  -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s)]])
        gamma_sz = gs*array([[(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r),  (1 - r)/2,  -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 - r)/2, -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r)]])

        # Reference 1, Equations 5.102
        B_gamma_MITC4 = zeros((2, 12))
        B_gamma_MITC4[0, :] = gamma_rz*sin(beta) - gamma_sz*sin(alpha)
        B_gamma_MITC4[1, :] = -gamma_rz*cos(beta) + gamma_sz*cos(alpha)

        # Return the [B] matrix for shear
        return B_gamma_MITC4

    #%%
    def B_m(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[s + 1, -s - 1, s - 1, -s + 1],
                                                  [r + 1, -r + 1, r - 1, -r - 1]]))

        # Reference 1, Example 5.5 (page 353)
        B_m = array([[dH[0, 0],    0,     dH[0, 1],    0,     dH[0, 2],    0,     dH[0, 3],    0    ],
                     [   0,     dH[1, 0],    0,     dH[1, 1],    0,     dH[1, 2],    0,     dH[1, 3]],
                     [dH[1, 0], dH[0, 0], dH[1, 1], dH[0, 1], dH[1, 2], dH[0, 2], dH[1, 3], dH[0, 3]]])

        return B_m

    #%%
    def Cb(self):
        '''
        Returns the stress-strain matrix for plate bending.
        '''

        # Referemce 1, Table 4.3, page 194
        nu = self.nu
        E = self.E
        h = self.t

        Cb = E*h**3/(12*(1 - nu**2))*array([[1,  nu,      0    ],
                                            [nu, 1,       0    ],
                                            [0,  0,  (1 - nu)/2]])

        return Cb

    #%%
    def Cs(self):
        '''
        Returns the stress-strain matrix for shear.
        '''
        # Reference 1, Equations (5.97), page 422
        k = 5/6
        E = self.E
        h = self.t
        nu = self.nu

        Cs = E*h*k/(2*(1 + nu))*array([[1, 0],
                                       [0, 1]])

        return Cs

    #%%
    def Cm(self):
        """
        Returns the stress-strain matrix for an isotropic or orthotropic plane stress element
        """

        # Apply the stiffness modification factors for each direction to obtain orthotropic
        # behavior. Stiffness modification factors of 1.0 in each direction (the default) will
        # model isotropic behavior. Orthotropic behavior is limited to the element's local
        # coordinate system.
        Ex = self.E*self.kx_mod
        Ey = self.E*self.ky_mod
        nu_xy = self.nu
        nu_yx = self.nu

        # The shear modulus will be unafected by orthotropic behavior
        # Logan, Appendix C.3, page 750
        G = self.E/(2*(1 + self.nu))

        # Gallagher, Equation 9.3, page 251
        Cm = 1/(1 - nu_xy*nu_yx)*array([[   Ex,    nu_yx*Ex,           0         ],
                                        [nu_xy*Ey,    Ey,              0         ],
                                        [    0,        0,     (1 - nu_xy*nu_yx)*G]])

        return Cm

    #%%
    def k_b(self):
        '''
        Returns the local stiffness matrix for bending stresses
        '''

        Cb = self.Cb()
        Cs = self.Cs()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the determinant of the Jacobian matrix for each gauss pointing
        # Doing this now will save us from doing it twice below
        J1 = det(self.J(gp, gp))
        J2 = det(self.J(-gp, gp))
        J3 = det(self.J(-gp, -gp))
        J4 = det(self.J(gp, -gp))

        # Get the bending B matrices for each gauss point
        B1 = self.B_kappa(gp, gp)
        B2 = self.B_kappa(-gp, gp)
        B3 = self.B_kappa(-gp, -gp)
        B4 = self.B_kappa(gp, -gp)

        # Create the stiffness matrix with bending stiffness terms
        # See Reference 1, Equation 5.94
        k_b = (matmul(B1.T, matmul(Cb, B1))*J1 +
             matmul(B2.T, matmul(Cb, B2))*J2 +
             matmul(B3.T, matmul(Cb, B3))*J3 +
             matmul(B4.T, matmul(Cb, B4))*J4)

        # Get the MITC4 shear B matrices for each gauss point
        B1 = self.B_gamma(gp, gp)
        B2 = self.B_gamma(-gp, gp)
        B3 = self.B_gamma(-gp, -gp)
        B4 = self.B_gamma(gp, -gp)

        # Alternatively the shear B matrix below could be used. However, this matrix is prone to
        # shear locking and will overestimate the stiffness.
        # B1 = self.B_gamma(gp, gp)
        # B2 = self.B_gamma(-gp, gp)
        # B3 = self.B_gamma(-gp, -gp)
        # B4 = self.B_gamma(gp, -gp)

        # # Add shear stiffness terms to the stiffness matrix
        k_s= (matmul(B1.T, matmul(Cs, B1))*J1 +
              matmul(B2.T, matmul(Cs, B2))*J2 +
              matmul(B3.T, matmul(Cs, B3))*J3 +
              matmul(B4.T, matmul(Cs, B4))*J4)

        k = k_b + k_s

        # Following Bathe's recommendation for the drilling degree of freedom
        # from Example 4.19 in "Finite Element Procedures, 2nd Ed.", calculate
        # the drilling stiffness as 1/1000 of the smallest diagonal term in
        # the element's stiffness matrix. This is not theoretically correct,
        # but it allows the model to solve without singularities, and should
        # have a minimal effect on the final solution. Bathe recommends 1/1000
        # as a value that is weak enough but not so small that it affect the
        # results. Bathe recommends looking at all the diagonals in the
        # combined bending plus membrane stiffness matrix. Some of those terms
        # relate to translational stiffness. It seems more rational to only
        # look at the terms relating to rotational stiffness. That will be
        # PyNite's approach.
        k_rz = min(abs(k[1, 1]), abs(k[2, 2]), abs(k[4, 4]), abs(k[5, 5]),
                   abs(k[7, 7]), abs(k[8, 8]), abs(k[10, 10]), abs(k[11, 11])
                   )/1000

        # Initialize the expanded stiffness matrix to all zeros
        k_exp = zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(12):

            # j = Unexpanded matrix column
            for j in range(12):

                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 3, 6, 9]:  # indices associated with deflection in z
                    m = 2*i + 2
                if i in [1, 4, 7, 10]:  # indices associated with rotation about x
                    m = 2*i + 1
                if i in [2, 5, 8, 11]:  # indices associated with rotation about y
                    m = 2*i

                # n = Expanded matrix column
                if j in [0, 3, 6, 9]:  # indices associated with deflection in z
                    n = 2*j + 2
                if j in [1, 4, 7, 10]:  # indices associated with rotation about x
                    n = 2*j + 1
                if j in [2, 5, 8, 11]:  # indices associated with rotation about y
                    n = 2*j

                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded
                # matrix
                k_exp[m, n] = k[i, j]

        # Add the drilling degree of freedom's weak spring
        k_exp[5, 5] = k_rz
        k_exp[11, 11] = k_rz
        k_exp[17, 17] = k_rz
        k_exp[23, 23] = k_rz

        return k_exp

    def area(self):
        gp = 1/3**0.5
        J1 = det(self.J(gp, gp))
        J2 = det(self.J(-gp, gp))
        J3 = det(self.J(-gp, -gp))
        J4 = det(self.J(gp, -gp))
        return J1 + J2 + J3 + J4

    #%%
    def k_m(self):
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        Cm = self.Cm()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the membrane B matrices for each gauss point
        # Doing this now will save us from doing it twice below
        B1 = self.B_m(gp, gp)
        B2 = self.B_m(-gp, gp)
        B3 = self.B_m(-gp, -gp)
        B4 = self.B_m(gp, -gp)

        # See reference 1 at the bottom of page 353, and reference 2 page 466
        k = t*(matmul(B1.T, matmul(Cm, B1))*det(self.J(gp, gp)) +
               matmul(B2.T, matmul(Cm, B2))*det(self.J(-gp, gp)) +
               matmul(B3.T, matmul(Cm, B3))*det(self.J(-gp, -gp)) +
               matmul(B4.T, matmul(Cm, B4))*det(self.J(gp, -gp)))

        k_exp = zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(8):

            # j = Unexpanded matrix column
            for j in range(8):

                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 2, 4, 6]:  # indices associated with displacement in x
                    m = i*3
                if i in [1, 3, 5, 7]:  # indices associated with displacement in y
                    m = i*3 - 2

                # n = Expanded matrix column
                if j in [0, 2, 4, 6]:  # indices associated with displacement in x
                    n = j*3
                if j in [1, 3, 5, 7]:  # indices associated with displacement in y
                    n = j*3 - 2

                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded matrix
                k_exp[m, n] = k[i, j]

        return k_exp

    def k(self):
        '''
        Returns the quad element's local stiffness matrix.
        '''
        self._local_coords()

        #node_positions = [
        #  np.array([0.0, 0.0,0.0]),
        #  np.array([1.5297058540778354, 0.0,0.0]),
        #  np.array([0.5491251783869153, 1.8434916702989301,0.0]),
        #  np.array([-0.6667948594698258, 0.7844645405527361,0.0]),
        #
        # print(self.x1, self.y1)
        # print(self.x2, self.y2)
        # print(self.x3, self.y3)
        # print(self.x4, self.y4)
        #
        # self.x1 = 0.0
        # self.y1 = 0.0
        # self.x2 = 1.5297058540778354
        # self.y2 = 0.0
        # self.x3 = 0.5491251783869153
        # self.y3 = 1.8434916702989301
        # self.x4 = -0.6667948594698258
        # self.y4 = 0.7844645405527361

        return self.k_b() + self.k_m()

        # # Recalculate the local coordinate system
        #
        # # Sum the bending and membrane stiffness matrices
        # return add(self.k_b(), self.k_m())

    def T(self):
        '''
        Returns the coordinate transformation matrix for the quad element.
        '''

        xi = self.i_node.X
        xj = self.j_node.X
        yi = self.i_node.Y
        yj = self.j_node.Y
        zi = self.i_node.Z
        zj = self.j_node.Z

        # Calculate the direction cosines for the local x-axis.The local x-axis will run from
        # the i-node to the j-node
        x = [xj - xi, yj - yi, zj - zi]

        # Divide the vector by its magnitude to produce a unit x-vector of
        # direction cosines
        mag = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
        x = [x[0]/mag, x[1]/mag, x[2]/mag]

        # The local y-axis will be in the plane of the plate. Find a vector in
        # the plate's local xy plane.
        xn = self.n_node.X
        yn = self.n_node.Y
        zn = self.n_node.Z
        xy = [xn - xi, yn - yi, zn - zi]

        # Find a vector perpendicular to the plate surface to get the
        # orientation of the local z-axis.
        z = cross(x, xy)

        # Divide the z-vector by its magnitude to produce a unit z-vector of
        # direction cosines.
        mag = (z[0]**2 + z[1]**2 + z[2]**2)**0.5
        z = [z[0]/mag, z[1]/mag, z[2]/mag]

        # Calculate the local y-axis as a vector perpendicular to the local z
        # and x-axes.
        y = cross(z, x)

        # Divide the y-vector by its magnitude to produce a unit vector of
        # direction cosines.
        mag = (y[0]**2 + y[1]**2 + y[2]**2)**0.5
        y = [y[0]/mag, y[1]/mag, y[2]/mag]

        # Create the direction cosines matrix.
        dirCos = array([x,
                        y,
                        z])

        # Build the transformation matrix.
        T = zeros((24, 24))
        T[0:3, 0:3] = dirCos
        T[3:6, 3:6] = dirCos
        T[6:9, 6:9] = dirCos
        T[9:12, 9:12] = dirCos
        T[12:15, 12:15] = dirCos
        T[15:18, 15:18] = dirCos
        T[18:21, 18:21] = dirCos
        T[21:24, 21:24] = dirCos

        # Return the transformation matrix.
        return T

    def K(self):
        '''
        Returns the quad element's global stiffness matrix
        '''

        # Get the transformation matrix
        T = self.T()

        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(inv(T), self.k()), T)



model = Model()
model.nodes[1] = Node(0.1,-0.5, 0.0 )
model.nodes[2] = Node(1.5, 0.3, 0.0 )
model.nodes[3] = Node(1.0, 1.2, 0.0 )
model.nodes[4] = Node(-0.5, 1.5, 0.0)
plate = MITC4('Plate 1', m_node=model.nodes[1],
                         n_node=model.nodes[2],
                         i_node=model.nodes[3],
                         j_node=model.nodes[4], t=3, material_name='Steel', model=model)
plate._local_coords()
import numpy as np

np.set_printoptions(precision=2,  linewidth=2000, suppress=True)
print(plate.K())