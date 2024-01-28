import numpy as np

class Solution:
    def __init__(self):
        self.loadcases = {}
    def get(self, loadcase, field):
        return self.loadcases[loadcase][field]

    def mises(self, stress):
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
        return np.sqrt(0.5 * ((sigma_x - sigma_y)**2 +
                              (sigma_y - sigma_z)**2 +
                              (sigma_z - sigma_x)**2 +
                              6 * (tau_xy**2 + tau_yz**2 + tau_zx**2)))

    def principal(self, stress):
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
        principal_stresses = np.zeros((sigma_x.shape[0], 3))
        for i in range(sigma_x.shape[0]):
            stress_tensor = np.array([[sigma_x[i], tau_xy[i], tau_zx[i]],
                                      [tau_xy[i], sigma_y[i], tau_yz[i]],
                                      [tau_zx[i], tau_yz[i], sigma_z[i]]])
            eigvals, _ = np.linalg.eig(stress_tensor)
            principal_stresses[i] = np.sort(eigvals)
        return principal_stresses

    def signed_mises(self, stress):
        mises_stress = self.mises(stress)
        principal_stresses_val = self.principal(stress)
        max_absolute_principal_stress = np.max(np.abs(principal_stresses_val), axis=1)
        sign = np.sign(np.sum(principal_stresses_val * (np.abs(principal_stresses_val) == max_absolute_principal_stress[:, None]), axis=1))
        return sign * mises_stress


    def list_fields(self, loadcase='1'):
        fields_dict = {}
        lc_fields = self.loadcases[str(loadcase)]
        for field, matrix in lc_fields.items():
            dim = matrix.shape[1]

            if dim == 6:
                if field == "STRESS":
                    fields_dict["mises"] = lambda f=field: self.mises(self.get(loadcase, f))
                    fields_dict["principal"] = lambda f=field: self.principal(self.get(loadcase, f))
                    fields_dict["signed_mises"] = lambda f=field: self.signed_mises(self.get(loadcase, f))

                # Check the field name and adjust suffixes accordingly
                if field == "DISPLACEMENT":
                    suffixes = ["_x", "_y", "_z"]
                else:
                    suffixes = ["_x", "_y", "_z", "_yz", "_zx", "_xy"]

                for idx, suffix in enumerate(suffixes):
                    fields_dict[field.lower() + suffix] = lambda f=field, i=idx: self.get(loadcase, f)[:,i]

                if field == "DISPLACEMENT":
                    fields_dict["displacement"] = lambda f=field: np.linalg.norm(self.get(loadcase, f)[:, :3], axis=1)
                    fields_dict["displacement_xyz"] = lambda f=field: self.get(loadcase, f)[:,:3]
                    fields_dict["rotation_x"] = lambda f=field: self.get(loadcase, f)[:,3]
                    fields_dict["rotation_y"] = lambda f=field: self.get(loadcase, f)[:,4]
                    fields_dict["rotation_z"] = lambda f=field: self.get(loadcase, f)[:,5]

            elif dim == 3:
                suffixes = ["_x", "_y", "_z"]
                for idx, suffix in enumerate(suffixes):
                    fields_dict[field.lower() + suffix] = lambda f=field, i=idx: self.get(loadcase, f)[:,i]
                fields_dict[field.lower()] = lambda f=field: np.linalg.norm(self.get(loadcase, f), axis=1)

            else:
                fields_dict[field.lower()] = lambda f=field: self.get(loadcase, f)

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
                    if current_loadcase:
                        loadcases[current_loadcase] = fields
                    fields = {}
                    current_loadcase = line.split()[1].strip()

                # Check for a field marker
                elif line.startswith("FIELD, NAME="):
                    if current_field_name:
                        fields[current_field_name] = np.array(matrix)
                    matrix = []
                    current_field_name = line.split("NAME=")[1].split(",")[0].strip()

                elif line == "END FIELD":
                    if current_field_name:
                        fields[current_field_name] = np.array(matrix)
                    current_field_name = None

                # If the line is a data row
                elif current_field_name:
                    temp = list(map(float, line.split()))
                    # replace nans with zeros
                    temp = [0 if np.isnan(x) else x for x in temp]
                    matrix.append(temp)

            if current_loadcase:
                loadcases[current_loadcase] = fields

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
