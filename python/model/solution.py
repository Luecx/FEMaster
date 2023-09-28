import numpy as np

class Solution:

    def __init__(self):
        self.loadcases = {}

    def mises(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
        mises = np.sqrt(0.5 * ((sigma_x - sigma_y) ** 2 + (sigma_y - sigma_z) ** 2 +
                               (sigma_z - sigma_x) ** 2 + 6 * (tau_xy ** 2 + tau_yz ** 2 + tau_zx ** 2)))
        return mises

    def principal(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
        num_points = sigma_x.shape[0]
        principal_stresses = np.zeros((num_points, 3))

        for i in range(num_points):
            stress_tensor = np.array([[sigma_x[i], tau_xy[i], tau_zx[i]],
                                      [tau_xy[i], sigma_y[i], tau_yz[i]],
                                      [tau_zx[i], tau_yz[i], sigma_z[i]]])
            eigvals, _ = np.linalg.eig(stress_tensor)
            principal_stresses[i] = np.sort(eigvals)

        return principal_stresses

    def signed_mises(self, loadcase='1'):
        mises_stress = self.mises(loadcase)
        principal_stresses_val = self.principal(loadcase)
        max_absolute_principal_stress = np.max(np.abs(principal_stresses_val), axis=1)
        sign = np.sign(np.sum(principal_stresses_val * (np.abs(principal_stresses_val) == max_absolute_principal_stress[:, None]), axis=1))
        signed_mises_stress = sign * mises_stress
        return signed_mises_stress

    def stress_x(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 0]

    def stress_y(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 1]

    def stress_z(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 2]

    def stress_yz(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 3]

    def stress_zx(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 4]

    def stress_xy(self, loadcase='1'):
        stress = self.loadcases[str(loadcase)]["STRESS"]
        return stress[:, 5]

    def displacement_x(self, loadcase='1'):
        displacement = self.loadcases[str(loadcase)]["DISPLACEMENT"]
        return displacement[:, 0]

    def displacement_y(self, loadcase='1'):
        displacement = self.loadcases[str(loadcase)]["DISPLACEMENT"]
        return displacement[:, 1]

    def displacement_z(self, loadcase='1'):
        displacement = self.loadcases[str(loadcase)]["DISPLACEMENT"]
        return displacement[:, 2]

    def displacement(self, loadcase='1'):
        displacement = self.loadcases[str(loadcase)]["DISPLACEMENT"]
        magnitude = np.linalg.norm(displacement, axis=1)
        return magnitude

    def list_fields(self, loadcase='1'):
        fields_dict = {}
        lc_fields = self.loadcases[str(loadcase)]
        for field, matrix in lc_fields.items():
            dim = matrix.shape[1]
            if dim == 6:
                if field == "STRESS":
                    fields_dict.update({
                        "mises": lambda: self.mises(loadcase),
                        "principal": lambda: self.principal(loadcase),
                        "signed_mises": lambda: self.signed_mises(loadcase),
                    })
                fields_dict.update({
                    field.lower() + "_x": lambda: self.get(loadcase, field)[:, 0],
                    field.lower() + "_y": lambda: self.get(loadcase, field)[:, 1],
                    field.lower() + "_z": lambda: self.get(loadcase, field)[:, 2],
                    field.lower() + "_yz": lambda: self.get(loadcase, field)[:, 3],
                    field.lower() + "_zx": lambda: self.get(loadcase, field)[:, 4],
                    field.lower() + "_xy": lambda: self.get(loadcase, field)[:, 5]
                })
            elif dim == 1:
                fields_dict[field.lower()] = lambda: self.get(loadcase, field)
            elif dim == 3:
                fields_dict[field.lower() + "_x"] = lambda: self.get(loadcase, field)[:, 0]
                fields_dict[field.lower() + "_y"] = lambda: self.get(loadcase, field)[:, 1]
                fields_dict[field.lower() + "_z"] = lambda: self.get(loadcase, field)[:, 2]
                fields_dict[field.lower()] = lambda: np.linalg.norm(self.get(loadcase, field), axis=1)
            else:
                pass  # handle other dimensions if needed

        return fields_dict

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
