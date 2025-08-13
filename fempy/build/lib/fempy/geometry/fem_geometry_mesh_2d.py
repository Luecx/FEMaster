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
    from .fem_geometry import Geometry
    import gmsh
    import numpy as np

    gmsh.initialize()
    gmsh.model.add("SegmentGroupMesh")
    gmsh.option.setNumber("General.Terminal", 0)

    curve_loops = []
    embedded_curves = []

    for group in segment_groups:
        boundary_points = group.get_points()
        if group.is_closed():
            boundary_points = boundary_points[:-1]

        print(f"boundary (closede={group.is_closed()})")
        for p in boundary_points:
            print(f"\t{p}")

        point_tags = []
        for pt in boundary_points:
            tag = gmsh.model.geo.addPoint(pt[0], pt[1], 0, tolerance)
            point_tags.append(tag)

        line_tags = []
        for i in range(len(point_tags) - 1):
            line = gmsh.model.geo.addLine(point_tags[i], point_tags[i + 1])
            line_tags.append(line)

        if group.is_closed():
            # Close the loop
            line = gmsh.model.geo.addLine(point_tags[-1], point_tags[0])
            line_tags.append(line)

            loop_tag = gmsh.model.geo.addCurveLoop(line_tags)
            curve_loops.append(loop_tag)
        else:
            embedded_curves.extend(line_tags)

    # Create surface
    surface_tag = gmsh.model.geo.addPlaneSurface(curve_loops)

    print("--------")
    print("Surface: ", surface_tag)

    gmsh.model.geo.synchronize()

    print("Synchronised")

    # Embed internal curves
    if embedded_curves:
        gmsh.model.mesh.embed(1, embedded_curves, 2, surface_tag)

    print("Embedded")

    # gmsh.option.setNumber("Mesh.Algorithm", 6)
    # gmsh.option.setNumber("Mesh.ElementOrder", 1 if not second_order else 2)
    #
    # if mesh_type in {1, 2}:
    #     gmsh.model.geo.mesh.setRecombine(2, surface_tag)

    gmsh.model.geo.synchronize()
    print("synchronised")
    gmsh.model.mesh.generate(2)
    print("generated")

    node_data = gmsh.model.mesh.getNodes()
    nodes = node_data[1].reshape((-1, 3))

    print(nodes)


    element_types, element_tags, node_tags_flattened = gmsh.model.mesh.getElements(dim=2)

    element_node_count = {
        2: 3,
        3: 4,
        9: 6,
        10: 9,
    }

    node_tags_per_element = []
    for elem_type, node_tags in zip(element_types, node_tags_flattened):
        num_nodes_per_elem = element_node_count.get(elem_type, 0)
        if num_nodes_per_elem > 0:
            node_tags_per_element.append(np.array(node_tags).reshape(-1, num_nodes_per_elem))

    geometry = Geometry(dimension=2)

    node_ids = {node_id: idx for idx, node_id in enumerate(node_data[0])}
    for idx, node in enumerate(nodes):
        geometry.add_node(idx + 1, node[0], node[1])

    elem_id = 1
    for elem_tags, nodes_per_elem in zip(element_tags, node_tags_per_element):
        for elem_tag, node_tags in zip(elem_tags, nodes_per_elem):
            node_indices = [node_ids[n] + 1 for n in node_tags]
            if len(node_tags) == 3:
                element_type = 'S3'
            elif len(node_tags) == 4:
                element_type = 'S4'
            elif len(node_tags) == 6:
                element_type = 'S6'
            elif len(node_tags) == 9:
                element_type = 'S8'
                node_indices = node_indices[:8]
            else:
                raise ValueError(f"Unsupported element type with {len(node_tags)} nodes.")

            geometry.add_element(elem_id, element_type, node_indices)
            elem_id += 1

    for group in segment_groups:
        geometry.add_node_set(group.name)
        for segment in group.segments:
            geometry.add_node_set(segment.name)

    for id, node in enumerate(geometry.nodes):
        if node is None:
            continue
        node_coords = np.asarray([node[0], node[1]])
        for group in segment_groups:
            for segment in group.segments:
                if segment.contains(node_coords):
                    geometry.add_node_to_set(segment.name, id)
                    geometry.add_node_to_set(group.name, id)

    gmsh.finalize()
    return geometry
