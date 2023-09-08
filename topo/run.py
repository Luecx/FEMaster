import subprocess
import os
from datetime import datetime
import numpy as np


def parse_output_file(filename):
    fields = {}
    current_field_name = None
    matrix = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            # Check for a field marker
            if line.startswith("FIELD, NAME="):
                # If we're currently processing a field, store its matrix
                if current_field_name:
                    fields[current_field_name] = np.array(matrix)

                # Reset the matrix and set the current field name
                matrix = []
                current_field_name = line.split("NAME=")[1].split(",")[0].strip()

            elif line == "END FIELD":
                if current_field_name:
                    fields[current_field_name] = np.array(matrix)
                current_field_name = None

            # If the line is a data row (and not a field marker or "END FIELD"), add it to the matrix
            elif current_field_name:
                matrix.append(list(map(float, line.split())))

    return fields


def run(model, densities):
    # Create the log directory if it doesn't exist
    if not os.path.exists(model.path + "/log"):
        os.makedirs(model.path + "/log")

    # Get the current timestamp and create the filename for the log file
    timestamp = datetime.now().strftime("%Y%m%d_%H-%M-%S")
    log_filename = f"{model.path}/log/solver_output_{timestamp}.log"

    # Write the input file
    model.write(densities, model.path + "/model.inp")

    # Run the solver and redirect both stdout and stderr to the log file
    with open(log_filename, 'w') as log_file:
        result = subprocess.run(["./solver", model.path + "/model.inp"], stdout=log_file, stderr=log_file)

    # Check for errors in execution
    if result.returncode != 0:
        with open(log_filename, 'r') as log_file:
            raise Exception("Error in executing the solver:", log_file.read())

    data = parse_output_file(model.path + '/model.inp.res')

    data['COMPLIANCE_ADJ'] = data['COMPLIANCE_ADJ']
    data['DENS_GRAD'] = data['DENS_GRAD']

    # os.remove(model.path + '/model.inp')
    # os.remove(model.path + '/model.inp.res')

    return data