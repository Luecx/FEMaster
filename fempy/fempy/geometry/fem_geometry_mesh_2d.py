
def mesh_interior(segment_groups, second_order=False, mesh_type=0, tolerance=1e-6):
    """
      Generate a 2D mesh for a given set of segment groups.

      Parameters
      ----------
      segment_groups : list of SegmentGroup
          List of SegmentGroup objects defining the boundaries of the mesh.
      second_order : bool, optional
          Whether to generate second-order elements. Default is False.
      mesh_type : int, optional
          Type of elements to generate:
              0 - Triangles only (default)
              1 - Mixed mesh (triangles and quadrilaterals)
              2 - Quadrilaterals only
      tolerance : float, optional
          Tolerance value for point comparison and boundary detection.

      Returns
      -------
      geometry : Geometry
          Geometry object containing the generated mesh nodes and elements.
    """
    # Lazy import to avoid circular dependency
    from .fem_geometry import Geometry
    import gmsh
    import numpy as np


    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("SegmentGroupMesh")
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress terminal output

    # Store the curve loops and surfaces
    curve_loops = []
    internal_curves = []

    # error if any of the segments are not closed
    # TODO: fix this
    for group in segment_groups:
        if not group.is_closed():
            raise ValueError(f"Segment group {group.name} is not closed")

    for group in segment_groups:
        # Generate boundary points from each segment group
        boundary_points = group.get_points()

        if group.is_closed():
            boundary_points = boundary_points[:-1]

        # Add points to Gmsh
        point_tags = []
        for pt in boundary_points:
            tag = gmsh.model.geo.addPoint(pt[0], pt[1], 0)
            point_tags.append(tag)

        if group.is_closed():
            # Add lines between consecutive points
            line_tags = []
            num_points = len(point_tags)
            for i in range(num_points):
                tag = gmsh.model.geo.addLine(point_tags[i], point_tags[(i + 1) % num_points])
                line_tags.append(tag)

            # Create a closed curve loop for this segment group
            curve_loop_tag = gmsh.model.geo.addCurveLoop(point_tags)
            curve_loops.append(curve_loop_tag)
        else:
            line_tags = []
            for i in range(len(point_tags) - 1):
                tag = gmsh.model.geo.addLine(point_tags[i], point_tags[(i + 1)])
                line_tags.append(tag)
            for i in range(len(line_tags) - 1):
                tag1 = line_tags[i]
                tag2 = line_tags[i + 1]
                gmsh.model.geo.addLine(tag1, tag2)

    surface_tag = gmsh.model.geo.addPlaneSurface(curve_loops)

    # Synchronize before meshing
    gmsh.model.geo.synchronize()

    # Set global meshing options
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Delaunay meshing for 2D
    gmsh.option.setNumber("Mesh.ElementOrder", 1 if not second_order else 2)  # First or second order elements

    # Handle mesh type: 0 = Triangles only, 1 = Mixed (tri + quad), 2 = Quad only
    if mesh_type in {1, 2}:  # Enable recombination for mixed or quad-only mesh
        gmsh.model.geo.mesh.setRecombine(2, surface_tag)

    # Synchronize again after recombination settings
    gmsh.model.geo.synchronize()

    # Generate the mesh
    gmsh.model.mesh.generate(2)

    # Extract nodes and elements
    node_data = gmsh.model.mesh.getNodes()
    nodes = node_data[1].reshape((-1, 3))

    element_types, element_tags, node_tags_flattened = gmsh.model.mesh.getElements(dim=2)

    # Manually define the number of nodes per element type
    element_node_count = {
        2: 3,  # 3-node triangle
        3: 4,  # 4-node quadrilateral
        9: 6,  # 6-node second-order triangle
        10: 9,  # 9-node second-order quadrilateral
    }

    # Convert flattened node tags to list of lists of node tags per element
    node_tags_per_element = []
    for elem_type, node_tags in zip(element_types, node_tags_flattened):
        num_nodes_per_elem = element_node_count.get(elem_type, 0)
        if num_nodes_per_elem > 0:
            node_tags_per_element.append(np.array(node_tags).reshape(-1, num_nodes_per_elem))

    # Initialize the resulting Geometry object
    geometry = Geometry(dimension=2)

    # Map node ids to indices
    node_ids = {node_id: idx for idx, node_id in enumerate(node_data[0])}

    # Add nodes to the geometry
    for idx, node in enumerate(nodes):
        geometry.add_node(idx + 1, node[0], node[1])

    elem_id = 1
    # Add elements to the geometry
    for elem_tags, nodes_per_elem in zip(element_tags, node_tags_per_element):
        for elem_tag, node_tags in zip(elem_tags, nodes_per_elem):
            node_indices = [node_ids[n] + 1 for n in node_tags]

            # Identify element type based on the number of nodes and create geometry
            if len(node_tags) == 3:
                element_type = 'S3'  # 3-node triangle
            elif len(node_tags) == 4:
                element_type = 'S4'  # 4-node quadrilateral
            elif len(node_tags) == 6:
                element_type = 'S6'  # 6-node second-order triangle
            elif len(node_tags) == 9:
                element_type = 'S8'  # Always convert 9-node to 8-node quadrilateral
                node_indices = node_indices[:8]
            else:
                raise ValueError(f"Unsupported element type with {len(node_tags)} nodes.")

            geometry.add_element(elem_id, element_type, node_indices)
            elem_id += 1

    # create a set for every group / segment in each group
    for group in segment_groups:
        geometry.add_node_set(group.name)
        for segment in group.segments:
            geometry.add_node_set(segment.name)

    # go through every node and check if it is in any segment
    for id, node in enumerate(geometry.nodes):
        if node is None:
            continue
        node_coords = np.asarray([node[0], node[1]])
        for group in segment_groups:
            for segment in group.segments:
                if segment.contains(node_coords):
                    geometry.add_node_to_set(segment.name, id)
                    geometry.add_node_to_set(group.name, id)


    # Finalize Gmsh
    gmsh.finalize()
    return geometry

