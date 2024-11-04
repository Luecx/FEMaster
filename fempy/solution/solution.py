import numpy as np

class Solution:
    def __init__(self):
        self.loadcases = {}

    def get(self, loadcase, field):
        return self.loadcases[loadcase][field]

    def list_fields(self, loadcase='1'):
        fields_dict = {}
        lc_fields = self.loadcases[loadcase]

        for field, matrix in lc_fields.items():
            dim = matrix.shape[1]

            # Defining suffixes
            suffixes = ["_x", "_y", "_z", "_yz", "_zx", "_xy"]

            # Create fields for available columns in the matrix
            for idx, suffix in enumerate(suffixes[:dim]):
                fields_dict[field.lower() + suffix] = lambda f=field, i=idx: self.get(loadcase, f)[:, i]

            # If 3 or more columns, add xyz and its magnitude
            if dim >= 3:
                fields_dict[field.lower() + "_xyz"] = lambda f=field: self.get(loadcase, f)[:, :3]
                fields_dict[field.lower() + "_xyz_mag"] = lambda f=field: np.linalg.norm(self.get(loadcase, f)[:, :3], axis=1)

            # If 6 or more columns, add yz, zx, xy and their magnitude
            if dim >= 6:
                fields_dict[field.lower() + "_yz_zx_xy"] = lambda f=field: self.get(loadcase, f)[:, 3:6]
                fields_dict[field.lower() + "_yz_zx_xy_mag"] = lambda f=field: np.linalg.norm(self.get(loadcase, f)[:, 3:6], axis=1)

            # Specific derived fields for STRESS
            if field == "STRESS":
                fields_dict["mises"] = lambda f=field: self.mises(self.get(loadcase, f))
                fields_dict["principal"] = lambda f=field: self.principal(self.get(loadcase, f))
                fields_dict["signed_mises"] = lambda f=field: self.signed_mises(self.get(loadcase, f))

            # Specific derived fields for DISPLACEMENT
            if field == "DISPLACEMENT":
                fields_dict["displacement"] = lambda f=field: np.linalg.norm(self.get(loadcase, f)[:, :3], axis=1)
                fields_dict["displacement_xyz"] = lambda f=field: self.get(loadcase, f)[:, :3]
                fields_dict["rotation_x"] = lambda f=field: self.get(loadcase, f)[:, 3]
                fields_dict["rotation_y"] = lambda f=field: self.get(loadcase, f)[:, 4]
                fields_dict["rotation_z"] = lambda f=field: self.get(loadcase, f)[:, 5]

        return fields_dict

    def list_fields_reduced(self, loadcase='1'):
        fields_dict = {}

        # Set the prefix for field names
        if loadcase is None:
            loadcases = self.loadcases.items()
        else:
            loadcases = [(loadcase, self.loadcases[loadcase])]

        for lc, lc_fields in loadcases:
            prefix = f"LC{lc}_" if loadcase is None else ""
            for field, matrix in lc_fields.items():
                # Specific derived fields for STRESS
                if field == "STRESS":
                    fields_dict[f"{prefix}mises"] = lambda f=field, lc=lc: self.mises(self.get(lc, f))
                    # fields_dict[f"{prefix}principal"] = lambda f=field, lc=lc: self.principal(self.get(lc, f))
                    # fields_dict[f"{prefix}signed_mises"] = lambda f=field, lc=lc: self.signed_mises(self.get(lc, f))

                    # fields_dict[f"{prefix}principal_dir_1"] = lambda f=field, lc=lc: self.principal_directions(self.get(lc, f))[:, 0,:]
                    # fields_dict[f"{prefix}principal_dir_2"] = lambda f=field, lc=lc: self.principal_directions(self.get(lc, f))[:, 1,:]
                    # fields_dict[f"{prefix}principal_dir_3"] = lambda f=field, lc=lc: self.principal_directions(self.get(lc, f))[:, 2,:]
                if field == "DISPLACEMENT":
                    fields_dict[f"{prefix}displacement_xyz"] = lambda f=field, lc=lc: self.get(lc, f)[:, :3]
                if "MODE_SHAPE" in field:
                    fields_dict[f"{prefix}{field.lower()}_xyz"] = lambda f=field, lc=lc: self.get(lc, f)[:, :3]
                # Store the raw data for other fields
                fields_dict[f"{prefix}{field.lower()}"] = lambda f=field, lc=lc: self.get(lc, f)

        return fields_dict


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

    def principal_directions(self, stress):
        """Compute the principal stress directions for each stress tensor."""
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T
        principal_stresses = np.zeros((sigma_x.shape[0], 3, 3))  # Shape: (num_elements, 3 principal directions, 3 components)

        for i in range(sigma_x.shape[0]):
            # Construct the stress tensor for the current element/point
            stress_tensor = np.array([[sigma_x[i], tau_xy[i], tau_zx[i]],
                                      [tau_xy[i], sigma_y[i], tau_yz[i]],
                                      [tau_zx[i], tau_yz[i], sigma_z[i]]])

            # Compute the eigenvalues and eigenvectors
            eigenvals, eigvecs = np.linalg.eig(stress_tensor)

            # Store the principal directions (eigenvectors)
            principal_stresses[i] = eigvecs * eigenvals

        return principal_stresses

    def signed_mises(self, stress):
        mises_stress = self.mises(stress)
        principal_stresses_val = self.principal(stress)
        max_absolute_principal_stress = np.max(np.abs(principal_stresses_val), axis=1)
        sign = np.sign(np.sum(principal_stresses_val * (np.abs(principal_stresses_val) == max_absolute_principal_stress[:, None]), axis=1))
        return sign * mises_stress

    @staticmethod
    def open(filename):
        import tqdm

        loadcases = {}
        current_loadcase = None
        current_field_name = None
        matrix = []

        with open(filename, 'r') as file:

            lines = file.readlines()
            # init progress bar
            pbar = tqdm.tqdm(total=len(lines), desc="Reading solution file")

            for line in lines:
                pbar.update(1)
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
                    # temp = [0 if np.isnan(x) else x for x in temp]
                    matrix.append(temp)

            if current_loadcase:
                loadcases[current_loadcase] = fields

        pbar.close()

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

# Example usage
# sol = Solution.open('path_to_your_file')
# print(sol)
# fields = sol.list_fields('1')
