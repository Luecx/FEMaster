
import sys
import os
import shutil
import subprocess
import numpy as np

from .filter import Filter, FilterFunction
from .adam import Adam
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

        # Run the solver and redirect both stdout and stderr to the log file
        with open(log_path, 'w') as log_file:
            result = subprocess.run([self.solver, model_path, "--ncpus", str(self.ncpus)], stdout=log_file, stderr=log_file)

        # if the residual is too large, exit
        with open(log_path, 'r') as log_file:
            for line in log_file:
                if "residual" in line.lower():
                    residual = float(line.split()[-1])
                    if residual > 1e-2:
                        raise Exception("Residual too large:", residual)

        # Check for errors in execution
        if result.returncode != 0:
            with open(log_path, 'r') as log_file:
                raise Exception("Error in executing the solver:", log_file.read())
        return self._read_output(iteration)

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

    def start(self, iterations):

        print("Starting optimisation")
        print("   elements in design set: ", np.sum(self.desi_mask))
        print("   target density        : ", self.target_density)
        print("   filter radius         : ", self.filter_radius)
        print("   symmetry radius       : ", self.symmetry_radius)
        print("   closest distance      : ", self.closest_distance)
        print("   move limit            : ", self.move_limit)
        print("   symmetry              : ", self.symmetries)

        # create density
        self.density = np.ones((len(self.desi_mask),))
        self.density[self.desi_mask] *= self.target_density

        # make sure all fields are initialised
        self.check_fields()

        # normal and symmetry filters
        grads_filter = Filter(sigma=self.filter_radius  , coords=self.mid_points,
                                     symmetries=self.symmetries, filter_func=FilterFunction.CONSTANT)
        value_filter = Filter(sigma=self.symmetry_radius, coords=self.mid_points,
                                     symmetries=self.symmetries, filter_func=FilterFunction.CONSTANT)

        # go through the iteration
        for iter in range(1, iterations+1):
            print("\n" + "-"*40)
            print("Iteration ", iter)
            print("-"*40)
            start_time = time.time()

            data = self._run(iter)
            # data = self._read_output(iter)
            run_time = time.time()
            vols = data[0]['VOLUME'].flatten()[self.desi_mask]
            dens = self.density               [self.desi_mask]

            # sensitivities merged from all loadcases
            comp    = np.sum([data[i]['COMPLIANCE_ADJ'  ].flatten()[self.desi_mask] for i in range(len(data))], axis=0)
            sens    = np.sum([data[i]['DENS_GRAD'       ].flatten()[self.desi_mask] for i in range(len(data))], axis=0)
            or_grad = np.sum([data[i]['ORIENTATION_GRAD']          [self.orientation_mask] for i in range(len(data))], axis=0)
            or_grad = -or_grad

            # store self.volumes
            self.volumes = vols
            sens = grads_filter.apply(sens)
            sens = np.minimum(0, sens)
            sens = sens / np.max(np.abs(sens))

            l1 = 1e-20
            l2 = 1e20

            while (l2-l1) / l1 > 1e-6:
                lm = (l2 + l1) / 2

                new_x = np.multiply(dens, np.sqrt(-sens/(lm)))
                new_x = np.minimum(dens + self.move_limit, new_x)
                new_x = np.minimum(1, new_x)
                new_x = np.maximum(dens - self.move_limit, new_x)
                new_x = np.maximum(self.min_density, new_x)
                new_x = self.custom_density_adjustment(new_x)

                if len(self.symmetries) > 0:
                    new_x = value_filter.apply(values=new_x)
                new_mass = np.multiply(new_x, self.volumes)
                mass_ratio = np.sum(new_mass) / np.sum(self.volumes)

                if np.abs(mass_ratio-self.target_density) < 1e-8:
                    break

                if mass_ratio - self.target_density > 0:
                    l1 = lm
                else:
                    l2 = lm

            change = np.linalg.norm(dens - new_x)
            self.density[self.desi_mask] = new_x

            # orientation optimization
            self.orientations[self.orientation_mask] = self.orientation_adam.step(or_grad)

            self.save(os.path.join(self.output_folder, 'iterations', str(iter), 'model.dat'))

            # self._clean_folder(iter)

            print("\nResults:")
            print(f"    Objective       : {np.sum(comp)}")
            print(f"    Change          : {change}")
            print(f"    Mass Constraint : {np.sum(np.multiply(new_x, vols)) / np.sum(vols)}")
            print(f"    Threshold       : {self.threshold()}")
            print(f"    Self time       : {time.time() - run_time}")
            print(f"    Solver time     : {run_time - start_time}")
            print(f"    Total time      : {time.time() - start_time}")
            print("Iteration ", iter, " completed.")

