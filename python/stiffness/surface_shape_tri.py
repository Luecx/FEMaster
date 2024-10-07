import sympy as sp

# Define symbols
r, s = sp.symbols('r s')

# Define shape functions
N1 = 1 - 3*(r + s) + 2*(r + s)**2
N2 = r*(2*r - 1)
N3 = s*(2*s - 1)
N4 = 4*r*(1 - r - s)
N5 = 4*r*s
N6 = 4*s*(1 - r - s)

# Calculate second derivatives
ddN1_r2 = sp.diff(N1, r, r)
ddN1_s2 = sp.diff(N1, s, s)
ddN1_rs = sp.diff(N1, r, s)

ddN2_r2 = sp.diff(N2, r, r)
ddN2_s2 = sp.diff(N2, s, s)
ddN2_rs = sp.diff(N2, r, s)

ddN3_r2 = sp.diff(N3, r, r)
ddN3_s2 = sp.diff(N3, s, s)
ddN3_rs = sp.diff(N3, r, s)

ddN4_r2 = sp.diff(N4, r, r)
ddN4_s2 = sp.diff(N4, s, s)
ddN4_rs = sp.diff(N4, r, s)

ddN5_r2 = sp.diff(N5, r, r)
ddN5_s2 = sp.diff(N5, s, s)
ddN5_rs = sp.diff(N5, r, s)

ddN6_r2 = sp.diff(N6, r, r)
ddN6_s2 = sp.diff(N6, s, s)
ddN6_rs = sp.diff(N6, r, s)

# Display results
print(f"Second derivatives for N1: ∂²N1/∂r² = {ddN1_r2}, ∂²N1/∂s² = {ddN1_s2}, ∂²N1/∂(rs) = {ddN1_rs}")
print(f"Second derivatives for N2: ∂²N2/∂r² = {ddN2_r2}, ∂²N2/∂s² = {ddN2_s2}, ∂²N2/∂(rs) = {ddN2_rs}")
print(f"Second derivatives for N3: ∂²N3/∂r² = {ddN3_r2}, ∂²N3/∂s² = {ddN3_s2}, ∂²N3/∂(rs) = {ddN3_rs}")
print(f"Second derivatives for N4: ∂²N4/∂r² = {ddN4_r2}, ∂²N4/∂s² = {ddN4_s2}, ∂²N4/∂(rs) = {ddN4_rs}")
print(f"Second derivatives for N5: ∂²N5/∂r² = {ddN5_r2}, ∂²N5/∂s² = {ddN5_s2}, ∂²N5/∂(rs) = {ddN5_rs}")
print(f"Second derivatives for N6: ∂²N6/∂r² = {ddN6_r2}, ∂²N6/∂s² = {ddN6_s2}, ∂²N6/∂(rs) = {ddN6_rs}")
