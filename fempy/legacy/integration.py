import sympy as sp
import pandas as pd

# Define variables
r, s, t = sp.symbols('r s t')

# Define the functions to integrate
functions = [
    sp.sympify(1.0),
    sp.sympify(2.0),
    sp.sympify(-1.0),
    sp.sympify(0.5),
    r,
    s,
    t,
    r + s + t,
    r * r - r * s - t * s + 3.0 * t * t,
    r * s + r * t + s * t,
    r * r + s * s + t * t,
    2.0 * r * r - 3.0 * s * s + 4.0 * t * s - r * s,
    r * r * r + s * s * s + t * t * t + r * s,
    r * r * r + s * s * t + t * t * r + s,
    r * r * r - 2.0 * s * s * s + 3.0 * r * s * t + r,
    r * r * s - s * t * t + 4.0 * r * s * t + t,
    r * r * r * r + s * s * s * s + t * t * t * t,
    r * r * r * r - s * s * s * s + 2.0 * r * r * s * s + t * t * t * t,
    r * r * s * s + s * s * t * t + t * t * r * r,
    r * r * r * r - 2.0 * r * r * s * s + s * s * s * s - 3.0 * r * s * t * t,
    r * r * r * r * r + s * s * s * s * s + t * t * t * t * t + r * r + s * s,
    r * r * r * s * s + s * s * s * t * t + t * t * t * r * r + r * r + s * s,
    r * r * r * r * r - 3.0 * r * r * r * s * s + t * t * r * s * s + r * r + s * s,
    r * r * s * s * s + s * s * s * t * t - t * t * t * r * r + r + s
]

# Defining the regions for different cases, with the specified order of integration
regions = {
    "line_a": [("r", -1, 1), ("s", 0), ("t", 0)],               # Case A
    "line_b": [("r", 0, 1), ("s", 0), ("t", 0)],                # Case B
    "triangle": [("r", 0, 1), ("s", 0, 1 - r), ("t", 0)],       # Case C
    "quad": [("r", -1, 1), ("s", -1, 1), ("t", 0)],             # Case D
    "hex": [("r", -1, 1), ("s", -1, 1), ("t", -1, 1)],          # Case E
    "tet": [("r", 0, 1), ("s", 0, 1 - r), ("t", 0, 1 - r - s)], # Case F
    "wedge": [("r", 0, 1), ("s", 0, 1 - r), ("t", -1, 1)],      # Case G
    "pyramid": [("t", 0, 1), ("r", -1+t, 1-t), ("s", -1+t, 1-t)]  # Case H
}

# Function to perform integration based on the specified order and constant variables
def integrate(function, limits):
    limits = limits[::-1]  # Reverse the limits to start from the innermost integral
    for limit in limits:
        var, *bounds = limit
        if len(bounds) == 2:  # If two bounds are provided, integrate
            function = sp.integrate(function, (var, bounds[0], bounds[1]))
        elif len(bounds) == 1:  # If one bound is provided, substitute as a constant
            function = function.subs(var, bounds[0])
    return function

# Recompute the table with the corrected approach
table = {}

for case_index in regions:
    region = regions[case_index]
    print(f"Case: {case_index}")
    table[case_index] = [integrate(f, region) for f in functions]

# Convert the results to a DataFrame and display
df = pd.DataFrame(table, index=[f"Function {i+1}" for i in range(len(functions))])

# Print the resulting DataFrame
print(df)
