import numpy as np
from .geometry import Geometry
from .solution import Solution
from scipy.optimize import curve_fit


def generate_beam(length, height, elem_length, elem_height, dent_func):
    geom = Geometry()
    n_elem_x = int(length // elem_length)
    n_elem_y = int(height // elem_height)

    geom.add_node_set('left_end')
    geom.add_node_set('right_end')
    geom.add_node_set('bottom_end')
    geom.add_node_set('top_end')
    geom.add_node_set('centerline')

    # Generate nodes
    node_id = 0
    for i in range(n_elem_x + 1):
        for j in range(n_elem_y + 1):
            x = i * elem_length
            y = j * elem_height - height / 2

            # Adjust y based on dent_func for the middle part
            factor = y / (height / 2)
            y += dent_func(x - length / 2) * factor
            geom.add_node(node_id, x, y)

            # Add to relevant sets
            if i == 0:
                geom.add_to_node_set('left_end', node_id)
            elif i == n_elem_x:
                geom.add_to_node_set('right_end', node_id)
            if j == 0:
                geom.add_to_node_set('bottom_end', node_id)
            elif j == n_elem_y:
                geom.add_to_node_set('top_end', node_id)
            if j == n_elem_y // 2:
                geom.add_to_node_set('centerline', node_id)

            node_id += 1

    # Generate C2D4 elements
    elem_id = 0
    for i in range(n_elem_x):
        for j in range(n_elem_y):
            bottom_left = i * (n_elem_y + 1) + j
            bottom_right = (i + 1) * (n_elem_y + 1) + j
            top_left = i * (n_elem_y + 1) + j + 1
            top_right = (i + 1) * (n_elem_y + 1) + j + 1

            # Add as C2D4
            node_ids = [bottom_left, bottom_right, top_right, top_left]
            geom.add_element(elem_id, 'C2D4', node_ids)
            elem_id += 1



    return geom

# Example usage:

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# Writing lines to the Input File
def append_lines(filename, lines):
    with open(filename, "a") as f:
        f.write("\n".join(lines))
        f.write("\n")

# Running the FEM Solver
def run_solver(input_filename):
    return subprocess.run(["./bin/FEMaster.exe", input_filename])

# Function to compute curvature radius based on dent function at x=0
def curvature_radius(D):
    f2_x0 = 2 / (D ** 2)
    R = 1 / np.abs(f2_x0)
    return R


# Main Code
if __name__ == '__main__':
    max_stresses = []
    R_values = []
    D_values = np.linspace(5, 0.1, 50)

    for D in D_values:
        filename = f"test_{D}.inp"

        # Dent function
        dent_function = lambda x: -1 * np.exp(-(x / D)**2)

        if D == 10:
            dent_function = lambda x: 0

        # Generate the geometry
        geom = generate_beam(20, 5, 0.3, 0.3, dent_func=dent_function)
        geom = geom.extrude(1, 1)
        geom.change_to_second_order()
        geom.write_input_deck(filename)

        # Append lines to the input deck
        lines_to_append = [
            "*SUPPORT, SUPPORT_COLLECTOR=SUPPS",
            "LEFT_END, 0, 0, 0",
            "*CLOAD, LOAD_COLLECTOR=LOADS",
            "RIGHT_END, 1, 0, 0",
            "*MATERIAL, NAME=MAT1",
            "*ELASTIC, TYPE=ISO",
            "210000, 0.3",
            "*SOLID SECTION, ELSET=EALL, MAT=MAT1",
            "*LOAD CASE, TYPE=LINEAR STATIC",
            "*LOAD",
            "LOADS",
            "*SUPPORT",
            "SUPPS",
            "*SOLVER, METHOD=DIRECT, DEVICE=CPU"
            "*END"
        ]

        append_lines(filename, lines_to_append)

        # Run the solver
        run_solver(filename)

        # Read the solution
        # Assuming you have a function Solution.open that returns the max Mises stress
        sol = Solution.open(f"test_{D}.inp.res")
        mises_stress = sol.mises(sol.get('1', 'STRESS'))
        max_stress = np.max(mises_stress)
        max_stresses.append(max_stress)

        # Compute the curvature radius R
        # Assuming a function curvature_radius(dent_function, D) that returns curvature radius for each D
        R = curvature_radius(D)

        R_values.append(R)

    # Curve Fitting
    def func(x, a, b):
        return a / (x ** b) + 1

    max_stresses = max_stresses[1:] / max_stresses[0]
    R_values = R_values[1:]

    params, _ = curve_fit(func, R_values, max_stresses)

    # Parameters: a
    a = params[0]
    b = params[1]

    # Generate fitted data for plotting
    R_values_fit = np.linspace(np.min(R_values), np.max(R_values), 400)
    max_stresses_fit = func(R_values_fit, a, b)

    # Plotting
    plt.figure()
    plt.scatter(R_values, max_stresses, label='Data', marker='o')
    plt.plot(R_values_fit, max_stresses_fit, label=f'Fit: a={a:.4f} b={b:.4f}' , linestyle='--')
    plt.xlabel("Curvature Radius R at x = 0")
    plt.ylabel("Max Mises Stress")
    plt.title("Max Mises Stress vs Curvature Radius")
    plt.grid(True)
    plt.legend()
    plt.show()