import numpy as np
import matplotlib.pyplot as plt
import time

def apply_symmetry(coords, values, symmetries):
    new_coords = coords.copy()
    new_values = values.copy()

    for symmetry, loc in symmetries.items():
        if symmetry == 'yz':
            mirrored_coords = new_coords.copy()
            mirrored_coords[:, 0] = 2*loc - mirrored_coords[:, 0]
            new_coords = np.vstack((new_coords, mirrored_coords))
            new_values = np.hstack((new_values, new_values))
        elif symmetry == 'xz':
            mirrored_coords = new_coords.copy()
            mirrored_coords[:, 1] = 2*loc - mirrored_coords[:, 1]
            new_coords = np.vstack((new_coords, mirrored_coords))
            new_values = np.hstack((new_values, new_values))
        elif symmetry == 'xy':
            mirrored_coords = new_coords.copy()
            mirrored_coords[:, 2] = 2*loc - mirrored_coords[:, 2]
            new_coords = np.vstack((new_coords, mirrored_coords))
            new_values = np.hstack((new_values, new_values))

    for symmetry, params in symmetries.items():
        if symmetry in ['x', 'y', 'z']:
            n_rotations, axis_loc = params
            for i in range(1, n_rotations):
                angle = 2 * np.pi * i / n_rotations
                rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                            [np.sin(angle), np.cos(angle)]])
                rotated_coords = new_coords.copy()
                if symmetry == 'x':
                    rotated_coords[:, 1:3] = (rotated_coords[:, 1:3] - axis_loc[1:]) @ rotation_matrix.T + axis_loc[1:]
                elif symmetry == 'y':
                    rotated_coords[:, ::2] = (rotated_coords[:, ::2] - axis_loc[::2]) @ rotation_matrix.T + axis_loc[::2]
                elif symmetry == 'z':
                    rotated_coords[:, :2] = (rotated_coords[:, :2] - axis_loc[:2]) @ rotation_matrix.T + axis_loc[:2]
                new_coords = np.vstack((new_coords, rotated_coords))
                new_values = np.hstack((new_values, new_values))

    return new_coords, new_values