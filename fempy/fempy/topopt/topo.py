
import sys
import os
import shutil
import subprocess
import numpy as np

from .filter import Filter, FilterFunction
from .adam import Adam
from .overhang import Overhang
from .plane import Plane
# from .viewer import Viewer
from ..geometry import Geometry
from ..solution import Solution
import time
import dill as pickle
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist


def compute_histogram(data, num_bins):
    # compute the histogram
    hist, bins = np.histogram(data, bins=num_bins)
    return hist, bins


class Optimiser:
    def __init__(self, input_deck, desi_set, loadcases, output_folder,
                 solver_path=None, method='indirect', device='gpu', exponent=2.5,
                 min_density=0.01, target_density=0.2, filter_radius=None,
                 symmetry_radius=None, move_limit=0.2, symmetries={}, custom_density_adjustment=None,
                 orientation_optimization=[False, False, False]):

        self.input_deck = input_deck
        self.desi_set = desi_set
        self.geometry = Geometry.read_input_deck(input_deck)
        self.volumes = np.zeros(len(self.geometry.elements))
        self.output_folder = output_folder

        self.loadcases = loadcases

        self.desi_mask = np.full(len(self.geometry.elements), False)
        self.desi_mask[self.geometry.elem_sets[self.desi_set]] = True

        mid_points_list = self.geometry.compute_element_midpoints()
        self.mid_points = np.array([mid_points_list[i] for i in range(len(mid_points_list)) if self.desi_mask[i] and mid_points_list[i] is not None])

        # solver settings
        self.solver = solver_path
        self.method = method
        self.device = device
        self.ncpus = 1

        # simp settings
        self.exponent = exponent
        self.min_density = min_density
        self.target_density = target_density
        self.filter_radius = filter_radius
        self.move_limit = move_limit

        # solution
        self.density = np.ones((len(self.desi_mask),))
        self.density[self.desi_mask] *= self.target_density

        # orientation optimization
        self.orientation_optimization = orientation_optimization
        self.orientation_mask = np.full((len(self.geometry.elements), 3), False)
        self.orientation_mask[self.desi_mask] = orientation_optimization
        self.orientations = np.zeros_like(self.orientation_mask, dtype=np.float64)
        self.orientations[self.orientation_mask] = np.random.uniform(-np.pi, np.pi, self.orientations[self.orientation_mask].shape)
        self.orientation_adam = Adam(learning_rate=self.move_limit, beta1=0.3, beta2=0.9, epsilon=1e-8, initial_values=
            self.orientations[self.orientation_mask])

        # constraints (symmetry)
        self.symmetries = symmetries
        self.symmetry_radius = symmetry_radius

        # fail if the element design set is not part of the geometry
        if self.desi_set not in self.geometry.elem_sets:
            raise ValueError("Design set is not part of geometry")

        # compute default radius in case its not specified
        min_dist = Filter(coords=self.mid_points, sigma=0).minimal_distance()
        self.closest_distance = min_dist
        if self.filter_radius is None:
            self.filter_radius = 10 * min_dist
        if self.symmetry_radius is None:
            self.symmetry_radius = 3 * min_dist

        # custom density adjustment
        if custom_density_adjustment is None:
           custom_density_adjustment = lambda x: x
        self.custom_density_adjustment = custom_density_adjustment

        # overhang constraint
        self.overhang_plane     = False
        self.overhang_dx        = False
        self.overhang_F         = False
        self.overhang_alpha     = False
        self.overhang_nodes     = False
        self.overhang_penalty   = False

    def check_fields(self):
        attrs = vars(self)
        none_fields = [attr for attr, value in attrs.items() if value is None]
        if none_fields:
            raise ValueError("The following fields are None: ", ", ".join(none_fields))

    def set_solver(self, path=None, method='indirect', device='gpu', ncpus=1):
        self.solver = path if path is not None else self.solver
        self.method = 'indirect' if method == 'indirect' else 'direct'
        self.device = 'gpu'      if device == 'gpu'      else 'cpu'
        self.ncpus  = ncpus

    def set_exponent(self, exponent=2.5):
        self.exponent = exponent

    def set_min_density(self, density=0.01):
        self.min_density = density

    def set_target_density(self, density=0.2):
        self.target_density = density

    def set_filter(self, radius):
        self.filter_radius = radius

    def set_move_limit(self, move_limit):
        self.move_limit = move_limit
        self.orientations_adam.learning_rate = move_limit

    def add_rotational_symmetry(self, axis, n_rotations, axis_loc):
        if axis in ['x', 'y', 'z']:
            self.symmetries[axis] = (n_rotations, axis_loc)
        else:
            raise ValueError("Invalid axis. Axis must be 'x', 'y', or 'z'.")

    def add_planar_symmetry(self, plane, loc):
        if plane in ['xy', 'yz', 'xz']:
            self.symmetries[plane] = loc
        else:
            raise ValueError("Invalid plane. Plane must be 'xy', 'yz', or 'xz'.")

    def set_symmetry_radius(self, radius):
        self.symmetry_radius = radius

    def _create_input(self, iteration):
        self.check_fields()

        # create output path
        output_path = os.path.join(self.output_folder, 'iterations', str(iteration))
        os.makedirs(output_path, exist_ok=True)

        # copy file
        new_file_path = os.path.join(output_path, 'model.inp')
        shutil.copy(self.input_deck, new_file_path)

        # append things to newly created file
        with open(new_file_path, 'a') as f:
            for lc in self.loadcases:
                load_cols = lc['load_cols']
                supp_cols = lc['supp_cols']
                f.write("\n")
                f.write("** AUTOMATICALLY GENERATED BY TOPOPT ** \n")
                f.write("*LOADCASE, TYPE=LINEAR STATIC TOPO\n")
                f.write("*LOAD\n")
                f.write(f"{','.join(load_cols)}\n")
                f.write("*SUPPORT\n")
                f.write(f"{','.join(supp_cols)}\n")
                f.write(f"*SOLVER, DEVICE={self.device}, METHOD={self.method}\n")
                f.write(f"*EXPONENT\n{self.exponent}\n")
                f.write("*DENSITY")
                for v in range(len(self.density)):
                    f.write(f"\n{v}, {self.density[v]}")
                f.write("\n*ORIENTATION")
                for v in range(len(self.orientations)):
                    f.write(f"\n{v}, {self.orientations[v][0]}, {self.orientations[v][1]}, {self.orientations[v][2]}")
                f.write("\n*END\n")

    def _run(self, iteration):
        self._create_input(iteration)
        output_path = os.path.join(self.output_folder, 'iterations', str(iteration))
        log_path = os.path.join(output_path, 'solver.log')
        model_path = os.path.join(output_path, 'model.inp')

        # Run the solver
        with open(log_path, 'w') as log_file:
            result = subprocess.run([self.solver, model_path, "--ncpus", str(self.ncpus)],
                                    stdout=log_file, stderr=log_file)

        # Check residual
        with open(log_path, 'r') as log_file:
            for line in log_file:
                if "residual" in line.lower():
                    residual = float(line.split()[-1])
                    if residual > 1e-2:
                        raise Exception("Residual too large:", residual)

        # Check exit code
        if result.returncode != 0:
            with open(log_path, 'r') as log_file:
                raise Exception("Error in executing the solver:", log_file.read())

        # Open .res once and return both raw and parsed results
        res_path = os.path.join(output_path, 'model.res')
        sol = Solution.open(res_path, loadingbar=False)
        results = [sol.loadcases[lc] for lc in sol.loadcases]
        return sol, results


    def _read_output(self, iteration):
        # create output path
        output_path = os.path.join(self.output_folder, 'iterations', str(iteration))

        # create file path
        file_path = os.path.join(output_path, 'model.res')

        # open file and return loadcase 1
        sol = Solution.open(file_path)
        results = []
        for lc in sol.loadcases:
            results.append(sol.loadcases[lc])
        return results

    def _clean_folder(self, iteration):
        # create output path
        output_path = os.path.join(self.output_folder, 'iterations', str(iteration))
        os.remove(os.path.join(output_path, 'model.inp'))
        os.remove(os.path.join(output_path, 'model.res'))

    def save(self, file_path):
        state = vars(self)
        with open(file_path, 'wb') as f:
            pickle.dump(state, f)

    def load(self, file_path):
        with open(file_path, 'rb') as f:
            state = pickle.load(f)
        for key, value in state.items():
            setattr(self, key, value)

    def load_it(self, iteration):
        file_path = os.path.join(self.output_folder, 'iterations', str(iteration), 'model.dat')
        self.load(file_path)

    def threshold(self):
        if self.volumes is None or np.all(self.volumes == 0):
            raise ValueError("Volumes have not been computed")
        # Define lower and upper bounds for the threshold
        lower_bound, upper_bound = 0, 1
        tol = 1e-6

        # Perform the bisection loop
        while upper_bound - lower_bound > tol:
            # Guess a threshold in between
            guess_threshold = (lower_bound + upper_bound) / 2

            # Count the number mass of the elements above
            mass_desi = np.sum(self.volumes)
            mass_dens = np.sum(self.volumes[(self.density[self.desi_mask] > guess_threshold)])

            # Adjust bounds based on the count
            if mass_dens / mass_desi > self.target_density:
                lower_bound = guess_threshold
            else:
                upper_bound = guess_threshold

        return (lower_bound + upper_bound) / 2

    def plot(self, threshold):
        v = Viewer()
        v.set_geometry(self.geometry)
        v.set_element_mask(self.density > threshold)
        v.set_data(type='element', data=self.density)
        v.set_data_range(0, 1)
        v.set_colorscheme('jet')
        v.coordinate_system()
        # v.set_boundaries()
        v.set_grid_xy()
        v.mainloop()

    def to_stl(self, threshold, file):
        v = Viewer()
        v.set_geometry(self.geometry)
        v.set_element_mask(self.density > threshold)
        v.write_to_stl(file)

    def set_overhang_constraint(self, plane, dx, F, alpha, maxnodes, penalty):
        self.overhang_plane = plane
        self.overhang_dx = dx
        self.overhang_F = F
        self.overhang_alpha = alpha
        self.overhang_nodes = maxnodes
        self.overhang_penalty = penalty

    def start(self, iterations):
        import time
        import os

        def _print_settings():
            print("Starting optimisation")
            print(f"   elements in design set: {np.sum(self.desi_mask)}")
            print(f"   target density        : {self.target_density}")
            print(f"   filter radius         : {self.filter_radius}")
            print(f"   symmetry radius       : {self.symmetry_radius}")
            print(f"   closest distance      : {self.closest_distance}")
            print(f"   move limit            : {self.move_limit}")
            print(f"   symmetry              : {self.symmetries}")

        def _initialize_density():
            self.density = np.ones(len(self.desi_mask))
            self.density[self.desi_mask] *= self.target_density
            self.check_fields()

        def _initialize_filters():
            grads_filter = Filter(self.mid_points, self.filter_radius, self.symmetries, FilterFunction.CONSTANT)
            value_filter = Filter(self.mid_points, self.symmetry_radius, self.symmetries, FilterFunction.CONSTANT)
            return grads_filter, value_filter

        def _initialize_overhang():
            if self.overhang_plane:
                return Overhang(
                    self.geometry, plane=self.overhang_plane,
                    F=self.overhang_F, alpha=self.overhang_alpha,
                    max_nodes=100000, element_set=self.desi_set
                )
            return None

        def _run_solver(iter):
            print("\n" + "-" * 40)
            print(f"Iteration {iter}")
            print("-" * 40)
            start_time = time.time()
            sol, data = self._run(iter)
            solver_end_time = time.time()
            return sol, data, start_time, solver_end_time

        def _compute_objectives(data):
            vols = data[0]["VOLUME"].flatten()[self.desi_mask]
            dens = self.density[self.desi_mask]
            comp = np.sum([d["COMPLIANCE_ADJ"].flatten()[self.desi_mask] for d in data], axis=0)
            sens = np.sum([d["DENS_GRAD"].flatten()[self.desi_mask] for d in data], axis=0)
            or_grad = -np.sum([d["ORIENTATION_GRAD"][self.orientation_mask] for d in data], axis=0)
            self.volumes = vols
            return comp, sens, or_grad, vols, dens

        def _apply_overhang_if_enabled(sens, comp_grad, overhang):
            if overhang:
                overhang_grad = overhang.loss_derivative(self.density, 1, self.overhang_penalty)[self.desi_mask]
                print(f"    Overhang gradient  : mean = {np.mean(overhang_grad):.6e}, var = {np.var(overhang_grad):.6e}")
                sens += overhang_grad
            else:
                overhang_grad = np.zeros_like(comp_grad)
            total_grad = sens.copy()
            return sens, overhang_grad, total_grad

        def _apply_sensitivities(sens, grads_filter):
            print(f"    Combined gradient  : mean = {np.mean(sens):.6e}, var = {np.var(sens):.6e}")
            sens = grads_filter.apply(sens)
            sens = np.minimum(0, sens)
            sens /= np.max(np.abs(sens))
            return sens

        def _update_densities(sens, dens, vols, value_filter):
            l1, l2 = 1e-20, 1e20
            while (l2 - l1) / l1 > 1e-6:
                lm = (l1 + l2) / 2
                new_x = np.multiply(dens, np.sqrt(-sens / lm))
                new_x = np.clip(new_x, dens - self.move_limit, dens + self.move_limit)
                new_x = np.clip(new_x, self.min_density, 1.0)
                new_x = self.custom_density_adjustment(new_x)
                if self.symmetries:
                    new_x = value_filter.apply(new_x)
                current_density = np.sum(new_x * vols) / np.sum(vols)
                if np.abs(current_density - self.target_density) < 1e-8:
                    break
                if current_density > self.target_density:
                    l1 = lm
                else:
                    l2 = lm
            return new_x

        def _finalize_density_update(new_x, dens):
            change = np.linalg.norm(dens - new_x)
            self.density[self.desi_mask] = new_x
            return change

        def _update_orientations(or_grad):
            self.orientations[self.orientation_mask] = self.orientation_adam.step(or_grad)

        def _write_results(sol, iter, comp_grad, overhang_grad, total_grad):
            field_len = len(self.desi_mask)
            comp_full = np.zeros(field_len)
            overhang_full = np.zeros(field_len)
            total_full = np.zeros(field_len)

            comp_full[self.desi_mask] = comp_grad
            overhang_full[self.desi_mask] = overhang_grad
            total_full[self.desi_mask] = total_grad

            sol.add_field("1", "COMP_GRAD", comp_full)
            sol.add_field("1", "CONSTRAINT_GRAD", overhang_full)
            sol.add_field("1", "TOTAL_GRAD", total_full)

            res_path = os.path.join(self.output_folder, 'iterations', str(iter), 'model.res')
            dat_path = os.path.join(self.output_folder, 'iterations', str(iter), 'model.dat')
            sol.write(res_path)
            self.save(dat_path)

        def _print_iteration_summary(comp, overhang, new_x, vols, change, start, solver_end):
            compliance_obj = np.sum(comp)
            if overhang:
                overhang_obj = overhang.loss(self.density, 1, self.overhang_penalty)
                total_obj = compliance_obj + overhang_obj
                print(f"\nResults:")
                print(f"    Objective       : {total_obj}")
                print(f"      Compliance    : {compliance_obj}")
                print(f"      Overhang      : {overhang_obj}")
            else:
                print(f"\nResults:")
                print(f"    Objective       : {compliance_obj}")
            print(f"    Change          : {change}")
            print(f"    Mass Constraint : {np.sum(new_x * vols) / np.sum(vols)}")
            print(f"    Threshold       : {self.threshold()}")
            print(f"    Self time       : {time.time() - solver_end:.2f}")
            print(f"    Solver time     : {solver_end - start:.2f}")
            print(f"    Total time      : {time.time() - start:.2f}")
            print(f"Iteration completed.\n")

        # ==== MAIN EXECUTION FLOW ====
        _print_settings()
        _initialize_density()
        grads_filter, value_filter = _initialize_filters()
        overhang = _initialize_overhang()

        for iter in range(1, iterations + 1):
            sol, data, start, solver_end = _run_solver(iter)
            comp, sens, or_grad, vols, dens = _compute_objectives(data)

            # Store raw compliance gradient before adding overhang
            comp_grad = sens.copy()
            sens, overhang_grad, total_grad = _apply_overhang_if_enabled(sens, comp_grad, overhang)

            # Filtering and clipping
            sens = _apply_sensitivities(sens, grads_filter)
            new_x = _update_densities(sens, dens, vols, value_filter)
            change = _finalize_density_update(new_x, dens)
            _update_orientations(or_grad)

            _write_results(sol, iter, comp_grad, overhang_grad, total_grad)
            _print_iteration_summary(comp, overhang, new_x, vols, change, start, solver_end)
