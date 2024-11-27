import copy
import numpy as np
from math import ceil

from ..geometry import Geometry

def generate_pipe(length, inner_radius, outer_radius, element_length, element_circum, element_radial, layers=None):
    geom = Geometry(dimension=2)
    n_elem_z = ceil(length // element_length)
    n_elem_u = ceil(2 * np.pi * inner_radius // element_circum)


    if layers is None:
        layers = {"LALL":(outer_radius - inner_radius)}
    n_elem_r = {key: ceil(layer // element_radial) for key, layer in layers.items()}

    # generate nodes.
    # generate inner layer first
    r_idx = 0
    for u in range(n_elem_u):
        id = r_idx * n_elem_u + u + 1
        x  = inner_radius * np.cos(u * 2 * np.pi / n_elem_u)
        y  = inner_radius * np.sin(u * 2 * np.pi / n_elem_u)
        geom.add_node(id, x, y, 0)
    r_idx += 1
    radius = inner_radius
    for layer in layers:
        r_start = radius
        r_end   = radius + layers[layer]

        for r in np.linspace(r_start, r_end, n_elem_r[layer] + 1)[1:]:
            for u in range(n_elem_u):
                id = r_idx * n_elem_u + u  + 1
                x = r * np.cos(u * 2 * np.pi / n_elem_u)
                y = r * np.sin(u * 2 * np.pi / n_elem_u)
                geom.add_node(id, x, y, 0)
            r_idx += 1
        radius = r_end
    # generate elements
    elem_id = 1
    r_idx = 0
    for layer, layer_r in layers.items():
        geom.add_element_set(layer)
        for r in range(n_elem_r[layer]):
            geom.add_element_set("LAYER_" + str(r_idx))
            for u in range(n_elem_u):
                node_ids = [r_idx * n_elem_u + u  + 1,
                            r_idx * n_elem_u + (u + 1) % n_elem_u  + 1,
                            (r_idx + 1) * n_elem_u + (u + 1) % n_elem_u + 1,
                            (r_idx + 1) * n_elem_u + u + 1]
                el_id = geom.add_element(elem_id, 'C2D4', node_ids)
                geom.add_element_to_set(layer, el_id)
                geom.add_element_to_set("LAYER_" + str(r_idx), el_id)
                elem_id += 1
            r_idx += 1

    geom = geom.extruded(n_elem_z, -length / n_elem_z)
    # geom = geom.as_second_order()

    return geom
    # generate the face and then extrude


    # n_elem_x = int(length // elements_length)
    # n_elem_y = int((outer_radius - inner_radius) // elements_radial)
    # n_elem_z = int(2 * np.pi * inner_radius // elements_circum)

geom = generate_pipe(length=500, element_length=1.5,
                     inner_radius=18, outer_radius=22,
                     element_circum=1, element_radial=0.5, layers={"L1": 1, "L2": 2, "L3": 1})
geom.write_input_deck("pipe.inp")