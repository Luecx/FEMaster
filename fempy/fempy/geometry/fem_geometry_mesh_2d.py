def mesh_interior(segment_groups, second_order=False, mesh_type=0, tolerance=1e-6):
    """
    Generate a 2D mesh for a given set of segment groups,
    allowing internal non-closed segments via OCC.fragment.

    Parameters
    ----------
    segment_groups : list of SegmentGroup
        List of SegmentGroup objects defining the boundaries
        and internal segments.
    second_order : bool, optional
        Whether to generate second-order elements. Default is False.
    mesh_type : int, optional
        Type of elements to generate:
            0 - Triangles only (default)
            1 - Mixed mesh (triangles and quads)
            2 - Quads only
    tolerance : float, optional
        Tolerance for point comparison (unused here but kept).

    Returns
    -------
    geometry : Geometry
        Geometry object containing mesh nodes and elements.
    """
    from .fem_geometry import Geometry
    import gmsh
    import numpy as np

    # ------------------------------------------------------------------
    # 1) Init Gmsh & OCC
    # ------------------------------------------------------------------
    gmsh.initialize()
    gmsh.model.add("SegmentGroupMesh")
    gmsh.option.setNumber("General.Terminal", 0)

    curve_loops     = []
    internal_curves = []

    # ------------------------------------------------------------------
    # 2) Build all curves (closed→loops, open→internal lines) in OCC
    # ------------------------------------------------------------------
    for group in segment_groups:
        pts = group.get_points()
        if group.is_closed():
            pts = pts[:-1]

        # add the points
        ptags = [gmsh.model.occ.addPoint(x, y, 0) for x, y in pts]

        # lines between consecutive points
        line_tags = []
        if group.is_closed():
            n = len(ptags)
            for i in range(n):
                j = (i + 1) % n
                line_tags.append(gmsh.model.occ.addLine(ptags[i], ptags[j]))
            # one OCC curve loop per closed group
            curve_loops.append(gmsh.model.occ.addCurveLoop(line_tags))
        else:
            # open chain → treat as “internal” cuts
            for i in range(len(ptags) - 1):
                tag = gmsh.model.occ.addLine(ptags[i], ptags[i + 1])
                internal_curves.append(tag)

    # ------------------------------------------------------------------
    # 3) Make the outer surface, then cut it by the internals
    # ------------------------------------------------------------------
    # single surface spanned by all closed loops
    surface_tag = gmsh.model.occ.addPlaneSurface(curve_loops)
    # push everything from OCC into the model
    gmsh.model.occ.synchronize()

    # fragment: cut that surface by each internal line
    gmsh.model.occ.fragment(
        [(2, surface_tag)],
        [(1, line_tag) for line_tag in internal_curves]
    )
    # push the newly created pieces back into the model
    gmsh.model.occ.synchronize()

    # ------------------------------------------------------------------
    # 4) Mesh options & generation
    # ------------------------------------------------------------------
    gmsh.option.setNumber("Mesh.Algorithm",    6)  # 2D Delaunay
    gmsh.option.setNumber("Mesh.ElementOrder", 1 if not second_order else 2)
    if mesh_type in {1, 2}:
        gmsh.option.setNumber("Mesh.RecombineAll", 1)     # recombine every surface

    # now actually generate the 2D mesh
    gmsh.model.mesh.generate(2)

    # ------------------------------------------------------------------
    # 5) Extract nodes & elements into your Geometry object
    # ------------------------------------------------------------------
    node_data = gmsh.model.mesh.getNodes()
    coords    = node_data[1].reshape((-1, 3))
    id_map    = {nid: i for i, nid in enumerate(node_data[0])}

    elt_types, elt_tags, elt_node_tags = gmsh.model.mesh.getElements(dim=2)
    node_counts = {2:3, 3:4, 9:6, 10:9}
    elems_per_type = []
    for etype, flat in zip(elt_types, elt_node_tags):
        n = node_counts.get(etype, 0)
        if n > 0:
            elems_per_type.append(np.array(flat).reshape(-1, n))

    geometry = Geometry(dimension=2)
    # add nodes
    for i, (x, y, _) in enumerate(coords):
        geometry.add_node(i+1, x, y)

    # add elements
    eid = 1
    for tags_list, nodes_arr in zip(elt_tags, elems_per_type):
        for _, nodelist in zip(tags_list, nodes_arr):
            inds = [id_map[n] + 1 for n in nodelist]
            cnt  = len(inds)
            if   cnt == 3: et = "S3"
            elif cnt == 4: et = "S4"
            elif cnt == 6: et = "S6"
            elif cnt == 9:
                et   = "S8"
                inds = inds[:8]
            else:
                raise ValueError(f"Unsupported element with {cnt} nodes")
            geometry.add_element(eid, et, inds)
            eid += 1

    # ------------------------------------------------------------------
    # 6) Build your node‐sets as before
    # ------------------------------------------------------------------
    for group in segment_groups:
        geometry.add_node_set(group.name)
        for seg in group.segments:
            geometry.add_node_set(seg.name)

    for idx, (x, y, _) in enumerate(coords):
        pt = np.array([x, y])
        for group in segment_groups:
            for seg in group.segments:
                if seg.contains(pt):
                    geometry.add_node_to_set(seg.name, idx)
                    geometry.add_node_to_set(group.name, idx)

    # ------------------------------------------------------------------
    # 7) Finish up
    # ------------------------------------------------------------------
    gmsh.finalize()
    return geometry
