from math import ceil

from fempy.solution.tovtk import Converter
from fempy.generate import generate_beam
from fempy.topopt import Optimiser

import numpy as np
import os

optimiser = Optimiser(
    input_deck="orient_topo.inp",
    desi_set="EALL",
    loadcases=[{'load_cols': ['loads'], 'supp_cols': ['supps']}],
    output_folder=f"./orient_topo/",
    solver_path="./bin/FEMaster",
    method='direct',
    device='cpu',
    exponent=2.5,
    min_density=0.01,
    target_density=0.25,
    filter_radius=3,
    symmetry_radius=0.5,
    move_limit=0.2,
    orientation_optimization=[False, False, True])


optimiser.start(10)