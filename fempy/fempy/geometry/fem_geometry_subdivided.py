import copy

def subdivided_geometry(geometry, n=1, only_quads=False):
    """
    Subdivide the given geometry into finer elements.
    
    Parameters
    ----------
    geometry : Geometry
        The Geometry object to be subdivided.
    n : int, optional
        Number of subdivision steps. Default is 1.
    only_quads : bool, optional
        Whether to convert triangles to quadrilaterals. Default is False.

    Returns
    -------
    Geometry
        Subdivided geometry object.
    """

    # Lazy import to avoid circular dependency
    from .fem_geometry import Geometry
    import numpy as np

    # Determine element order
    element_order = geometry.determine_element_order()
    if element_order == 0:
        raise ValueError("Element order could not be determined or no elements are present.")
    if element_order == 2:
        raise ValueError("Subdivision is not supported for second-order elements.")

    # Create a copy of the current geometry
    result_geom = Geometry(geometry.dimension)
    result_geom.nodes = geometry.nodes.copy()
    result_geom.node_sets = {name: ids.copy() for name, ids in geometry.node_sets.items()}
    result_geom.elem_sets = {name: ids.copy() for name, ids in geometry.elem_sets.items()}
    result_geom.elements = copy.deepcopy(geometry.elements)

    # Apply n subdivision steps
    for _ in range(n):
        next_geom = Geometry(result_geom.dimension)
        next_geom.nodes = result_geom.nodes.copy()
        next_geom.node_sets = {name: ids.copy() for name, ids in result_geom.node_sets.items()}
        next_geom.elem_sets = {name: ids.copy() for name, ids in result_geom.elem_sets.items()}
        next_geom.elements = copy.deepcopy(result_geom.elements)

        midpoints = next_geom.compute_edge_midpoints()

        # Subdivide each element
        for elem in next_geom.elements:
            if elem is not None:
                elem.subdivide(midpoints, next_geom, only_quads=only_quads)

        result_geom = next_geom

    return result_geom
