import sympy as sp
import numpy as np

# Define the symbols
n1, n2, n3 = sp.symbols('n1 n2 n3')

# Define the shape functions
phi1 = 1/4 * (n1 + n2 - 1) * ((1 + n1) * (1 + n2) - n3 + (n1 * n2 * n3) / (1 - n3))
phi2 = 1/4 * (-n1 + n2 - 1) * ((1 - n1) * (1 + n2) - n3 - (n1 * n2 * n3) / (1 - n3))
phi3 = 1/4 * (-n1 - n2 - 1) * ((1 - n1) * (1 - n2) - n3 + (n1 * n2 * n3) / (1 - n3))
phi4 = 1/4 * (n1 - n2 - 1) * ((1 + n1) * (1 - n2) - n3 - (n1 * n2 * n3) / (1 - n3))
phi5 = n3 * (2 * n3 - 1)
phi6 = ((1 + n1 - n3) * (1 - n1 - n3) * (1 + n2 - n3)) / (2 * (1 - n3))
phi7 = ((1 + n2 - n3) * (1 - n2 - n3) * (1 - n1 - n3)) / (2 * (1 - n3))
phi8 = ((1 + n1 - n3) * (1 - n1 - n3) * (1 - n2 - n3)) / (2 * (1 - n3))
phi9 = ((1 + n2 - n3) * (1 - n2 - n3) * (1 + n1 - n3)) / (2 * (1 - n3))
phi10 = n3 * (1 + n1 - n3) * (1 + n2 - n3) / (1 - n3)
phi11 = n3 * (1 - n1 - n3) * (1 + n2 - n3) / (1 - n3)
phi12 = n3 * (1 - n1 - n3) * (1 - n2 - n3) / (1 - n3)
phi13 = n3 * (1 + n1 - n3) * (1 - n2 - n3) / (1 - n3)

# List of shape functions
phis = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10, phi11, phi12, phi13]

# Define the node coordinates
nodes = [
    (1, 1, 0),   # Node 1
    (-1, 1, 0),  # Node 2
    (-1, -1, 0), # Node 3
    (1, -1, 0),  # Node 4
    (0, 0, 1),   # Node 5
    (0, 1, 0),  # Node 6 (midpoint between nodes 1 and 2)
    (-1, 0, 0), # Node 7 (midpoint between nodes 2 and 3)
    (0, -1, 0), # Node 8 (midpoint between nodes 3 and 4)
    (1, 0, 0), # Node 9 (midpoint between nodes 4 and 1)
    (1/2, 1/2, 1/2), # Node 10 (midpoint between nodes 1 and 5)
    (-1/2, 1/2, 1/2), # Node 11 (midpoint between nodes 2 and 5)
    (-1/2, -1/2, 1/2), # Node 12 (midpoint between nodes 3 and 5)
    (1/2, -1/2, 1/2) # Node 13 (midpoint between nodes 4 and 5)
]

# Create a matrix to store evaluations of phi at each node
evaluation_matrix = sp.zeros(len(nodes), len(phis))

print(len(nodes))
print(len(phis))

# Evaluate each phi at each node and store in the matrix
for i, node in enumerate(nodes):
    n1_val, n2_val, n3_val = node
    for j, phi in enumerate(phis):
        n1_val_adjusted = n1_val + 1e-10
        n2_val_adjusted = n2_val + 1e-10
        n3_val_adjusted = n3_val + 1e-10
        evaluation_matrix[i, j] = phi.subs({n1: n1_val_adjusted, n2: n2_val_adjusted, n3: n3_val_adjusted})

# Compute the derivatives of each shape function with respect to n1, n2, and n3
derivatives = []
for phi in phis:
    d_phi_n1 = sp.diff(phi, n1)
    d_phi_n2 = sp.diff(phi, n2)
    d_phi_n3 = sp.diff(phi, n3)
    derivatives.append((d_phi_n1, d_phi_n2, d_phi_n3))

# Print the derivatives
for i, (d_phi_n1, d_phi_n2, d_phi_n3) in enumerate(derivatives):
    print(f"Derivatives of phi{i+1}:")
    print(f"  d_phi/dn1: {d_phi_n1}")
    print(f"  d_phi/dn2: {d_phi_n2}")
    print(f"  d_phi/dn3: {d_phi_n3}")
    print()

sp.pprint(evaluation_matrix)
