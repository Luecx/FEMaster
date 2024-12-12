from sympy import symbols, Matrix, diff, Function, expand
import sympy as sp

# Define the variables
r, s, t = symbols('r s t')
gr = Matrix(symbols('grx gry grz'))
gs = Matrix(symbols('gsx gsy gsz'))
gt = Matrix(symbols('gtx gty gtz'))

# Define the u, h, and v vectors (position specific)
u = [Matrix(symbols(f'u_{i}x u_{i}y u_{i}z')) for i in range(1, 5)]
h = [symbols(f'h_{i}') for i in range(1, 5)]
v = [Matrix(symbols(f'v_{i}x v_{i}y v_{i}z')) for i in range(1, 5)]
theta = [Matrix(symbols(f'theta_{i}x theta_{i}y theta_{i}z')) for i in range(1, 5)]

# Define the N functions
N = [Function(f'N{i}')(r, s) for i in range(1, 5)]


def cross(vec1, vec2):
    return Matrix([
        vec1[1] * vec2[2] - vec1[2] * vec2[1],
        vec1[2] * vec2[0] - vec1[0] * vec2[2],
        vec1[0] * vec2[1] - vec1[1] * vec2[0]
    ])

htv1 = h[0] * cross(theta[0], v[0])
htv2 = h[1] * cross(theta[1], v[1])
htv3 = h[2] * cross(theta[2], v[2])
htv4 = h[3] * cross(theta[3], v[3])

# concat horizontally
u_mat   = Matrix.hstack(u[0], u[1], u[2], u[3])
htv_mat = Matrix.hstack(htv1, htv2, htv3, htv4)

N_mat   = Matrix([N[0], N[1], N[2], N[3]])
gr_mat  = Matrix([gr])
gs_mat  = Matrix([gs])
gt_mat  = Matrix([gt])

u = (u_mat + t / 2 * htv_mat) * N_mat

e_11 = gr_mat.dot(u.diff(r))
e_22 = gs_mat.dot(u.diff(s))
e_12 = gr_mat.dot(u.diff(s)) + gs_mat.dot(u.diff(r))
e_23 = gt_mat.dot(u.diff(s)) + gs_mat.dot(u.diff(t))
e_13 = gt_mat.dot(u.diff(r)) + gr_mat.dot(u.diff(t))

u_vars = [f"u_{i+1}{comp}" for i in range(1) for comp in ['x', 'y', 'z']]
theta_vars = [f"theta_{i+1}{comp}" for i in range(1) for comp in ['x', 'y', 'z']]
my_vars = u_vars + theta_vars


print("---- e11 ----")
for var in my_vars:
    coeffs = e_11.diff(symbols(var))
    print(var, coeffs)
print ("---- e22 ----")
for var in my_vars:
    coeffs = e_22.diff(symbols(var))
    print(var, coeffs)

print ("---- e12 ----")
for var in my_vars:
    coeffs = e_12.diff(symbols(var))
    print(var, coeffs)
print("---- e23 ----")
for var in my_vars:
    coeffs = e_23.diff(symbols(var))
    print(var, coeffs)

print("---- e13 ----")
for var in my_vars:
    coeffs = e_13.diff(symbols(var))
    print(var, coeffs)

# Fully expand the expression for u

# e11_expanded = expand(e_11)
# print(e11_expanded)
#
# u_terms = {}
# theta_terms = {}
#
# for i in range(4):
#     for comp in ['x', 'y', 'z']:
#         u_key = f"u_{i+1}{comp}"
#         theta_key = f"theta_{i+1}{comp}"
#         u_terms[u_key] = e11_expanded.collect(symbols(u_key))
#         theta_terms[theta_key] = e11_expanded.collect(symbols(theta_key))
#
#
# for key, val in u_terms.items():
#     print(f"{key}: {val}")