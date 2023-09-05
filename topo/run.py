
import subprocess
import os
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
    # Write the input file
    model.write(densities, "model.inp")

    print(densities[0:15])

    # Run the solver
    result = subprocess.run(["./solver", "model.inp"], capture_output=True, text=True)

    # Check for errors in execution
    if result.returncode != 0:
        raise Exception("Error in executing the solver:", result.stderr)

    data = parse_output_file('model.inp.res')

    data['COMPLIANCE_ADJ'] = data['COMPLIANCE_ADJ'] * 100
    data['DENS_GRAD'] = data['DENS_GRAD'] * 100
    data['DENS_GRAD'] = model.filter(data['DENS_GRAD'] )

    print(sum(data['COMPLIANCE_ADJ']))

    os.remove('model.inp')
    os.remove('model.inp.res')

    return data