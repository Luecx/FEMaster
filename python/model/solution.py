import numpy as np

class Solution:

    def __init__(self):
        self.loadcases = {}

    @staticmethod
    def open(filename):
        loadcases = {}
        current_loadcase = None
        current_field_name = None
        matrix = []

        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()

                # Check for a loadcase marker
                if line.startswith("LC"):
                    # If we're currently processing a loadcase, store its fields
                    if current_loadcase:
                        loadcases[current_loadcase] = fields

                    # Reset the fields and set the current loadcase
                    fields = {}
                    current_loadcase = line.split()[1].strip()

                # Check for a field marker
                elif line.startswith("FIELD, NAME="):
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

                # If the line is a data row (and not a field marker, "END FIELD", or loadcase marker), add it to the matrix
                elif current_field_name:
                    matrix.append(list(map(float, line.split())))

            # Store the last loadcase or field
            if current_loadcase:
                loadcases[current_loadcase] = fields
            elif current_field_name:
                fields[current_field_name] = np.array(matrix)

        sol = Solution()
        sol.loadcases = loadcases
        return sol

    def __str__(self):
        output = ""
        for loadcase, fields in self.loadcases.items():
            output += f"LC: {loadcase}\n"
            for field_name, field_matrix in fields.items():
                output += f"  - Field: {field_name}, Dimensions: {field_matrix.shape}\n"
        return output
