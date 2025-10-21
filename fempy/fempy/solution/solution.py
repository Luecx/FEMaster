import re
import numpy as np

# ---------- Helpers for series detection & naming ----------

# Erkenne eine numerische Endung (optional mit _ oder - davor)
# Beispiele:  "U0", "U_1", "disp-u-23", "stress_002"
_NUM_SUFFIX = re.compile(r"^(?P<base>.*?)(?:[_\-]?)(?P<idx>\d+)$", re.IGNORECASE)

def split_series(name: str):
    """
    Zerlegt Feldnamen in (basis, index), falls am Ende eine Zahl steht.
    Rückgabe: (basis_string_ohne_index, int_index oder None)
    """
    m = _NUM_SUFFIX.match(name)
    if not m:
        return name, None
    base = m.group("base")
    try:
        idx = int(m.group("idx"))
    except ValueError:
        idx = None
    return base, idx

def join_series(base: str, idx):
    """Fügt Basis + _ + idx wieder zusammen, falls idx != None."""
    return f"{base}_{idx}" if idx is not None else base

def is_stress_name(name: str) -> bool:
    return "stress" in name.lower()

def looks_like_displacement_name(name: str) -> bool:
    n = name.lower()
    return ("displacement" in n) or ("_u_" in n) or (n.endswith("_u")) or (n == "u")

def replace_first_u_token_with_xyz(base: str) -> str:
    """
    Ersetze das erste Vorkommen von '_u_' oder '-u-' durch '_u_xyz_' bzw. '-u_xyz-'.
    Falls keines gefunden wird, hänge zur Not '_xyz' an.
    """
    if "_u_" in base:
        return base.replace("_u_", "_u_xyz_", 1)
    if "-u-" in base:
        return base.replace("-u-", "-u_xyz-", 1)
    if base.lower().endswith("_u"):
        return base + "_xyz"
    return base + "_xyz"

def _collapse_indices(indices):
    """
    [0,1,2,5,7,8,9] -> '0…2, 5, 7…9'
    """
    if not indices:
        return ""
    indices = sorted(set(indices))
    ranges = []
    start = prev = indices[0]
    for x in indices[1:]:
        if x == prev + 1:
            prev = x
            continue
        # schließe aktuellen block
        ranges.append((start, prev))
        start = prev = x
    ranges.append((start, prev))
    parts = []
    for a, b in ranges:
        parts.append(str(a) if a == b else f"{a}…{b}")
    return ", ".join(parts)

def _collapse_series_keys(keys):
    """
    Nimmt eine Liste von Namen (z. B. ['U_0','U_1','U_2','V_7']) und
    liefert kollabierte Darstellung:
        {'U': '0…2', 'V': '7'}
    """
    groups = {}
    for k in keys:
        base, idx = split_series(k)
        if idx is None:
            groups.setdefault(base, set())  # leere Menge → '—'
        else:
            groups.setdefault(base, set()).add(idx)
    collapsed = {}
    for base, idxs in groups.items():
        if not idxs:
            collapsed[base] = "—"
        else:
            collapsed[base] = _collapse_indices(sorted(idxs))
    return collapsed


class Solution:
    def __init__(self):
        self.loadcases = {}

    def get(self, loadcase, field):
        return self.loadcases[loadcase][field]

    # --------------------------------------------------------------------------
    # Volle Feldliste (ohne Reduktion), behält alte API bei.
    # Erzeugt komponenten-skalare (_x/_y/_z …), xyz-Tripel, Magnituden etc.
    # --------------------------------------------------------------------------
    def list_fields(self, loadcase='1'):
        fields_dict = {}
        lc_fields = self.loadcases[loadcase]

        for field, matrix in lc_fields.items():
            dim = matrix.shape[1]

            # Komponenten-Suffixe (für 6er Spannungen)
            suffixes = ["_x", "_y", "_z", "_yz", "_zx", "_xy"]

            # Einzelkomponenten (so weit verfügbar)
            for idx, suffix in enumerate(suffixes[:dim]):
                fields_dict[field.lower() + suffix] = (lambda f=field, i=idx:
                                                       self.get(loadcase, f)[:, i])

            # Ab 3 Spalten -> xyz + Betrag
            if dim >= 3:
                fields_dict[field.lower() + "_xyz"] = (lambda f=field:
                                                       self.get(loadcase, f)[:, :3])
                fields_dict[field.lower() + "_xyz_mag"] = (lambda f=field:
                                                           np.linalg.norm(self.get(loadcase, f)[:, :3], axis=1))

            # Ab 6 Spalten -> Schubtripel + Betrag
            if dim >= 6:
                fields_dict[field.lower() + "_yz_zx_xy"] = (lambda f=field:
                                                            self.get(loadcase, f)[:, 3:6])
                fields_dict[field.lower() + "_yz_zx_xy_mag"] = (lambda f=field:
                                                                np.linalg.norm(self.get(loadcase, f)[:, 3:6], axis=1))

            # Spezielle Ableitungen für STRESS
            if field.upper() == "STRESS":
                fields_dict["mises"] = (lambda f=field: self.mises(self.get(loadcase, f)))
                fields_dict["principal"] = (lambda f=field: self.principal(self.get(loadcase, f)))
                fields_dict["signed_mises"] = (lambda f=field: self.signed_mises(self.get(loadcase, f)))

            # Spezielles für DISPLACEMENT
            if field.upper() == "DISPLACEMENT":
                fields_dict["displacement"] = (lambda f=field:
                                               np.linalg.norm(self.get(loadcase, f)[:, :3], axis=1))
                fields_dict["displacement_xyz"] = (lambda f=field:
                                                   self.get(loadcase, f)[:, :3])
                fields_dict["rotation_x"] = (lambda f=field: self.get(loadcase, f)[:, 3])
                fields_dict["rotation_y"] = (lambda f=field: self.get(loadcase, f)[:, 4])
                fields_dict["rotation_z"] = (lambda f=field: self.get(loadcase, f)[:, 5])

        return fields_dict

    # --------------------------------------------------------------------------
    # Reduzierte Felder für Exporter (zeitserien-freundlich):
    # - erkennt Endziffern als Zeitindex und erzeugt konsistente Namen
    # - für *stress*-Felder (6 Spalten): <base>_mises[_idx]
    # - für "displacement"/"_u_" etc.: <base>_xyz[_idx] (nur die ersten 3 Komponenten)
    # - für Mode/Buckling: *_xyz[_idx]
    # - rohes Feld immer zusätzlich als <prefix><field.lower()>[_idx]
    # --------------------------------------------------------------------------
    def list_fields_reduced(self, loadcase='1'):
        fields_dict = {}

        # Mehrere LCs? Dann voranstellen: LC<id>_
        if loadcase is None:
            loadcases = self.loadcases.items()
        else:
            loadcases = [(loadcase, self.loadcases[loadcase])]

        for lc, lc_fields in loadcases:
            lc_prefix = f"LC{lc}_" if loadcase is None else ""

            for orig_name, matrix in lc_fields.items():
                base_no_idx, idx = split_series(orig_name)            # Basis ohne Zeitindex
                out_suffix = f"_{idx}" if idx is not None else ""     # _<idx> oder ""

                # Normalisierte Basen für abgeleitete Felder
                base_lower = base_no_idx.lower()
                name_raw   = f"{lc_prefix}{base_lower}{out_suffix}"   # Rohdaten-Name

                # 1) Rohdaten immer anbieten
                fields_dict[name_raw] = (lambda f=orig_name, l=lc: self.get(l, f))

                # 2) STRESS: von Mises, wenn 6 Komponenten
                if is_stress_name(base_no_idx) and matrix.ndim == 2 and matrix.shape[1] >= 6:
                    mises_name = f"{lc_prefix}{base_lower}_mises{out_suffix}"
                    fields_dict[mises_name] = (lambda f=orig_name, l=lc: self.mises(self.get(l, f)))

                # 3) DISPLACEMENT / *_u_* / U: nur erste 3 Komponenten als Vektor
                if matrix.ndim == 2 and matrix.shape[1] >= 3 and looks_like_displacement_name(base_no_idx):
                    # Name wie gewünscht: z.B. "..._u_xyz_<idx>"
                    xyz_base = replace_first_u_token_with_xyz(base_no_idx)
                    xyz_name = f"{lc_prefix}{xyz_base.lower()}{out_suffix}"
                    fields_dict[xyz_name] = (lambda f=orig_name, l=lc: self.get(l, f)[:, :3])

                    # dazu noch eine „kanonische“ Displacement-Vec (hilfreich für ParaView Warp)
                    disp_name = f"{lc_prefix}displacement{out_suffix}"
                    fields_dict[disp_name] = (lambda f=orig_name, l=lc: self.get(l, f)[:, :3])

                    # Magnitude optional nützlich
                    mag_name = f"{lc_prefix}{xyz_base.lower()}_mag{out_suffix}"
                    fields_dict[mag_name] = (lambda f=orig_name, l=lc:
                                             np.linalg.norm(self.get(l, f)[:, :3], axis=1))

                # 4) Mode-Shapes & Buckling-Shapes: *_xyz
                if matrix.ndim == 2 and matrix.shape[1] >= 3 and (
                        "mode_shape" in base_lower or "buckling_mode" in base_lower):
                    mode_xyz = f"{lc_prefix}{base_lower}_xyz{out_suffix}"
                    fields_dict[mode_xyz] = (lambda f=orig_name, l=lc: self.get(l, f)[:, :3])

        return fields_dict

    # ------------------------ Derived measures for stress ------------------------

    def mises(self, stress):
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T[:6]
        return np.sqrt(0.5 * ((sigma_x - sigma_y)**2 +
                              (sigma_y - sigma_z)**2 +
                              (sigma_z - sigma_x)**2 +
                              6 * (tau_xy**2 + tau_yz**2 + tau_zx**2)))

    def principal(self, stress):
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T[:6]
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
        sigma_x, sigma_y, sigma_z, tau_yz, tau_zx, tau_xy = stress.T[:6]
        principal_stresses = np.zeros((sigma_x.shape[0], 3, 3))  # (n, dir, comp)
        for i in range(sigma_x.shape[0]):
            stress_tensor = np.array([[sigma_x[i], tau_xy[i], tau_zx[i]],
                                      [tau_xy[i], sigma_y[i], tau_yz[i]],
                                      [tau_zx[i], tau_yz[i], sigma_z[i]]])
            eigenvals, eigvecs = np.linalg.eig(stress_tensor)
            principal_stresses[i] = eigvecs * eigenvals
        return principal_stresses

    def signed_mises(self, stress):
        mises_stress = self.mises(stress)
        principal_stresses_val = self.principal(stress)
        max_abs = np.max(np.abs(principal_stresses_val), axis=1)
        sign = np.sign(np.sum(principal_stresses_val * (np.abs(principal_stresses_val) == max_abs[:, None]), axis=1))
        return sign * mises_stress

    # ------------------------ IO ------------------------

    def add_field(self, loadcase, field_name, data):
        """
        Adds or overrides a field in the specified loadcase.

        Parameters:
        -----------
        loadcase : str or int
            Loadcase name or number.
        field_name : str
            Name of the field (e.g. 'SUPPORT' or 'DISPLACEMENT_u_7', 'STRESS_12').
            Eine Endziffer wird später als Zeitindex erkannt.
        data : np.ndarray
            Array shape (n,) oder (n,d).
        """
        loadcase = str(loadcase)
        if loadcase not in self.loadcases:
            self.loadcases[loadcase] = {}
        self.loadcases[loadcase][field_name] = np.atleast_2d(data).T if data.ndim == 1 else data

    def write(self, filepath):
        """
        Writes the current Solution to a result file in FEMaster format.
        """
        with open(filepath, 'w') as f:
            for lc, fields in self.loadcases.items():
                f.write(f"LC {lc}\n")
                for name, data in fields.items():
                    f.write(f"FIELD, NAME={name}\n")
                    for row in data:
                        f.write(" ".join(f"{v:.6e}" for v in row) + "\n")
                    f.write("END FIELD\n")

    @staticmethod
    def open(filename, loadingbar=True):
        import tqdm

        loadcases = {}
        current_loadcase = None
        current_field_name = None
        matrix = []

        with open(filename, 'r') as file:
            lines = file.readlines()
            if loadingbar:
                pbar = tqdm.tqdm(total=len(lines), desc="Reading solution file")

            for line in lines:
                if loadingbar:
                    pbar.update(1)
                line = line.strip()

                if line.startswith("LC"):
                    if current_loadcase:
                        loadcases[current_loadcase] = fields
                    fields = {}
                    current_loadcase = line.split()[1].strip()

                elif line.startswith("FIELD, NAME="):
                    if current_field_name:
                        fields[current_field_name] = np.array(matrix)
                    matrix = []
                    current_field_name = line.split("NAME=")[1].split(",")[0].strip()

                elif line == "END FIELD":
                    if current_field_name:
                        fields[current_field_name] = np.array(matrix)
                    current_field_name = None

                elif current_field_name:
                    temp = list(map(float, line.split()))
                    matrix.append(temp)

            if current_loadcase:
                loadcases[current_loadcase] = fields

        if loadingbar:
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

    # ------------------------ Reporting / Summary ------------------------

    def summary(self, include_reduced: bool = True) -> str:
        """
        Erzeugt einen Textbericht über alle Loadcases, Originalfelder (+ Shapes),
        und – optional – abgeleitete/reduzierte Felder. Serien werden kollabiert
        angezeigt (z. B. 'U_xyz_[0…100]').

        Parameters
        ----------
        include_reduced : bool
            Wenn True, werden zusätzlich alle Namen aus list_fields_reduced()
            ausgewertet und gruppiert (ohne Arrays zu berechnen).

        Returns
        -------
        str : schön formatierter Bericht.
        """
        lines = []
        lines.append("=== Solution Summary ===")
        lcs = list(self.loadcases.keys())
        lines.append(f"Loadcases: {', '.join(lcs) if lcs else '(none)'}")
        lines.append("")

        # 1) Pro LC: Originalfelder + Shapes
        for lc, fields in self.loadcases.items():
            lines.append(f"[LC {lc}] Original fields:")
            if not fields:
                lines.append("  (none)")
            else:
                # Gruppiere auch die Originalnamen nach Basis/Serie
                orig_by_base = {}
                for name, arr in fields.items():
                    base, idx = split_series(name)
                    shape = tuple(arr.shape)
                    entry = orig_by_base.setdefault(base, {"idxs": [], "shapes": set()})
                    if idx is not None:
                        entry["idxs"].append(idx)
                    entry["shapes"].add(shape)

                for base, info in sorted(orig_by_base.items()):
                    idxs = _collapse_indices(info["idxs"])
                    shapes = ", ".join(sorted(map(str, info["shapes"])))
                    if idxs:
                        lines.append(f"  - {base}[{idxs}]  shapes: {shapes}")
                    else:
                        lines.append(f"  - {base}          shapes: {shapes}")

            lines.append("")

            # 2) Reduzierte/abgeleitete Felder (nur Namen, Serien kollabiert)
            if include_reduced:
                reduced = self.list_fields_reduced(lc)
                # nur Namen abgreifen
                keys = list(reduced.keys())
                groups = _collapse_series_keys(keys)
                if groups:
                    lines.append(f"[LC {lc}] Reduced / derived fields (collapsed):")
                    for base, rng in sorted(groups.items()):
                        suffix = f"[{rng}]" if rng != "—" else ""
                        lines.append(f"  - {base}{suffix}")
                else:
                    lines.append(f"[LC {lc}] Reduced / derived fields: (none)")
                lines.append("")

        return "\n".join(lines)
