import os
import numpy as np
import matplotlib.pyplot as plt
from fempy.solution.solution import Solution

# Define the template for the input deck
INPUT_TEMPLATE = """
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
*ELEMENT, TYPE=C3D8
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=mat1
*ELASTIC, TYPE=ORTHOTROPIC
76262, 6655, 6655 , 2429, 4374, 4374, 0.37 , 0.33, 0.33

*MATERIAL, NAME=mat2
*ELASTIC, TYPE=ISOTROPIC
100, 0.3

*SOLID SECTION, ELSET=EALL, MATERIAL=mat1

*SUPPORT, SUPPORT_COLLECTOR=SUPPS
1, 0, 0, 0
4, 0,  , 0
5, 0,  ,
8, 0,  ,
*CLOAD, LOAD_COLLECTOR=LOADS
2, 1, 0, 0
3, 1, 0, 0
6, 1, 0, 0
7, 1, 0, 0

*LOAD CASE, type = linear static topo
*load
LOADS
*support
SUPPS
*DENSITY
1, 1
*ORIENTATION
1, 0, 0, {angle}
*END
"""

# Path to FEM solver
FEM_SOLVER_PATH = "./bin/FEMaster"

# Directory to store generated input files and results
OUTPUT_DIR = "input_decks"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def generate_input_deck(angle, filename):
    """
    Generates an input deck file with the specified orientation angle.
    """
    with open(filename, "w") as f:
        f.write(INPUT_TEMPLATE.format(angle=angle))

def call_fem_solver(input_file):
    """
    Calls the FEM solver with the specified input file, suppressing its output.
    """
    with open(os.devnull, 'w') as devnull:
        os.system(f"{FEM_SOLVER_PATH} {input_file} > /dev/null 2>&1")

def read_results(result_file):
    """
    Reads a .res file using the Solution class and extracts compliance and gradient data.
    """
    sol = Solution.open(result_file)
    angle_grad = sol.list_fields()['angle_grad_xyz']()
    compliance = sol.list_fields()['compliance_raw_x']()

    elem_id = 1  # Assuming a single element with ID 1
    return angle_grad[elem_id], compliance[elem_id]

def main():
    # Initialize lists to store results
    angles = []
    compliances = []
    gradients = []

    # Generate angles from 0 to 2π in steps of 0.1
    angle_values = np.arange(0, 2 * np.pi + 0.1, 0.1)

    for i, angle in enumerate(angle_values):
        angle = round(angle, 2)  # Round to avoid floating-point quirks
        input_filename = os.path.join(OUTPUT_DIR, f"input_{angle:.2f}.inp")
        result_filename = os.path.join(OUTPUT_DIR, f"input_{angle:.2f}.res")

        # Generate input deck
        generate_input_deck(angle, input_filename)

        # Print the current index
        print(f"Processing angle index: {i}, angle: {angle:.2f}")

        # Run FEM solver
        call_fem_solver(input_filename)

        # Read and store results
        if os.path.exists(result_filename):
            angle_grad, compliance = read_results(result_filename)
            angles.append(angle)
            compliances.append(compliance)
            gradients.append(angle_grad)

        # Clean up input and result files
        os.remove(input_filename)
        if os.path.exists(result_filename):
            os.remove(result_filename)

    # Convert results to numpy arrays for easier manipulation
    angles = np.array(angles)
    compliances = np.array(compliances)
    gradients = np.array(gradients)

    # Plot compliance vs. angles on the left y-axis and normalized gradients on the right y-axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Compliance on the left y-axis
    ax1.plot(angles, compliances, marker="o", label="Compliance", color="blue")
    ax1.set_xlabel("Angle (radians)")
    ax1.set_ylabel("Compliance", color="blue")
    ax1.tick_params(axis='y', labelcolor="blue")

    # Orientation gradient on the right y-axis
    ax2 = ax1.twinx()
    ax2.plot(angles, gradients, marker="x", label="Normalized Gradient", color="orange")
    ax2.set_ylabel("Orientation Gradient", color="orange")
    ax2.tick_params(axis='y', labelcolor="orange")

    # Title and grid
    plt.title("Compliance and Normalized Orientation Gradient vs. Angle")
    ax1.grid(True)

    # Show plot
    plt.show()

if __name__ == "__main__":
    main()
