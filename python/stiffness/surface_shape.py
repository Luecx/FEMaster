import sympy as sp

# Define symbolic variables
r, s = sp.symbols('r s')

# Define the shape functions
N1 = 0.25 * (1 - r) * (1 - s) * (-1 - r - s)
N2 = 0.25 * (1 + r) * (1 - s) * (-1 + r - s)
N3 = 0.25 * (1 + r) * (1 + s) * (-1 + r + s)
N4 = 0.25 * (1 - r) * (1 + s) * (-1 - r + s)
N5 = 0.5 * (1 - r*r) * (1 - s)
N6 = 0.5 * (1 + r) * (1 - s*s)
N7 = 0.5 * (1 - r*r) * (1 + s)
N8 = 0.5 * (1 - r) * (1 - s*s)

# Store shape functions in a matrix
N = sp.Matrix([N1, N2, N3, N4, N5, N6, N7, N8])

# Compute first derivatives
dN_dr = N.diff(r)
dN_ds = N.diff(s)

# Compute second derivatives
d2N_drr = dN_dr.diff(r)
d2N_dss = dN_ds.diff(s)
d2N_drs = dN_dr.diff(s)

# Display the results
print("Shape Functions (N):")
sp.pprint(sp.simplify(N))
print("\nFirst Derivatives w.r.t r (dN/dr):")
sp.pprint(sp.simplify(dN_dr))
print("\nFirst Derivatives w.r.t s (dN/ds):")
sp.pprint(sp.simplify(dN_ds))
print("\nSecond Derivatives w.r.t r (d2N/drr):")
sp.pprint(sp.simplify(d2N_drr))
print("\nSecond Derivatives w.r.t s (d2N/dss):")
sp.pprint(sp.simplify(d2N_dss))
print("\nSecond Derivatives w.r.t rs (d2N/drs):")
sp.pprint(sp.simplify(d2N_drs))
