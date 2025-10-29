from __future__ import annotations
"""
model_mesh.py — gmsh-backed 2D meshing for SegmentGroups

Builds native Gmsh geometry (Line, CircleArc, Spline, BSpline) from your
segments, hands subdivision to Gmsh (Transfinite), generates a 2D mesh,
and maps it back into your Model (nodes, elements, nodesets).

Special:
- If the outer closed group is a rectangle of 4 StraightSegments with right
  angles and no holes, we create a **perfect quad grid** via Transfinite Surface.

Element mapping:
- 3-node tri   -> C2D3
- 4-node quad  -> C2D4
- 6-node tri   -> C2D6
- 9-node quad  -> C2D8 (center node dropped → serendipity 8-node)
"""

from typing import Sequence, Tuple, Dict, List, TYPE_CHECKING, Optional
import numpy as np

# type-only to avoid circular import at import time
if TYPE_CHECKING:
    from .model import Model

from ..nodes.node import Node
from ..sets.nodeset import NodeSet
from ..elements.elements import Elements
from ..geometry import (
    StraightSegment,
    CircleSegment,
    CurvedSegment,
    BSplineSegment,
    Filet,
)

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def _rounded_key(x: float, y: float, eps: float = 1e-12) -> Tuple[float, float]:
    """Snap (x, y) to an epsilon grid to get a stable dict key (robust endpoint sharing)."""
    return (round(float(x) / eps) * eps, round(float(y) / eps) * eps)

def _add_point_cached(occ, cache: Dict[Tuple[float, float], int], p) -> int:
    """Create or reuse an OCC Point for 2D coordinate p (p can be (x,y) or np.ndarray)."""
    px = float(p[0]); py = float(p[1])
    key = _rounded_key(px, py)
    tag = cache.get(key)
    if tag is None:
        tag = int(occ.addPoint(px, py, 0.0))
        cache[key] = tag
    return tag

def _circle_center_from_three(p1, p2, p3) -> Tuple[float, float]:
    """Compute circle center from three 2D points."""
    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    x3, y3 = float(p3[0]), float(p3[1])
    a = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2
    b = (x1**2 + y1**2) * (y3 - y2) + (x2**2 + y2**2) * (y1 - y3) + (x3**2 + y3**2) * (y2 - y1)
    c = (x1**2 + y1**2) * (x2 - x3) + (x2**2 + y2**2) * (x3 - x1) + (x3**2 + y3**2) * (x1 - x2)
    cx = -b / (2.0 * a)
    cy = -c / (2.0 * a)
    return float(cx), float(cy)

# -----------------------------------------------------------------------------
# Curves per segment (collect transfinite wishes; don't set yet)
# -----------------------------------------------------------------------------

def _add_curve_for_segment(occ, pcache: Dict[Tuple[float, float], int], seg) -> Tuple[int, Optional[int]]:
    """Create a native OCC curve for seg and return (curve_tag, n_subdivisions_or_None)."""
    nsub = getattr(seg, "n_subdivisions", None)
    if isinstance(seg, StraightSegment):
        a = _add_point_cached(occ, pcache, seg.p_start)
        b = _add_point_cached(occ, pcache, seg.p_end)
        return int(occ.addLine(a, b)), nsub
    if isinstance(seg, CircleSegment):
        a = _add_point_cached(occ, pcache, seg.p_start)
        b = _add_point_cached(occ, pcache, seg.p_end)
        c = _add_point_cached(occ, pcache, seg.center)
        return int(occ.addCircleArc(a, c, b)), nsub
    if isinstance(seg, Filet):
        a = _add_point_cached(occ, pcache, seg.p_start)
        b = _add_point_cached(occ, pcache, seg.p_end)
        mid = seg.at(0.5)
        cx, cy = _circle_center_from_three(seg.p_start, mid, seg.p_end)
        c = _add_point_cached(occ, pcache, (cx, cy))
        return int(occ.addCircleArc(a, c, b)), nsub
    if isinstance(seg, CurvedSegment):
        a = _add_point_cached(occ, pcache, seg.p_start)
        m = _add_point_cached(occ, pcache, np.asarray(seg.mid))
        b = _add_point_cached(occ, pcache, seg.p_end)
        return int(occ.addSpline([a, m, b])), nsub
    if isinstance(seg, BSplineSegment):
        p_tags = [_add_point_cached(occ, pcache, np.asarray(p)) for p in seg.control_points]
        if getattr(seg, "closed", False) and len(p_tags) > 0:
            p_tags = list(p_tags) + [p_tags[0]]  # robust closure
        return int(occ.addBSpline(list(map(int, p_tags)))), nsub

    # Fallback: polyline
    pts = [_add_point_cached(occ, pcache, p) for p in seg.get_points()]
    if len(pts) < 2:
        raise ValueError("Segment produced < 2 points; cannot create curve.")
    first = int(occ.addLine(int(pts[0]), int(pts[1])))
    for s, e in zip(pts[1:-1], pts[2:]):
        occ.addLine(int(s), int(e))
    return int(first), nsub

# -----------------------------------------------------------------------------
# Rectangle detection for perfect quads
# -----------------------------------------------------------------------------

def _is_rectangular_group(
        group, tol_len: float = 1e-9, tol_orth: float = 1e-9
) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[List[int]], Optional[List[int]]]:
    """Detect if a closed group is a rectangle of 4 StraightSegments.

    Returns (is_rect, corner_coords, horiz_idx, vert_idx).
    corner_coords are start points of segments in group order.
    """
    segs = list(group.segments)
    if len(segs) != 4:
        return False, None, None, None
    for s in segs:
        if not isinstance(s, StraightSegment):
            return False, None, None, None

    corners = [(float(s.p_start[0]), float(s.p_start[1])) for s in segs]
    end0 = (float(segs[-1].p_end[0]), float(segs[-1].p_end[1]))
    if abs(corners[0][0] - end0[0]) > 1e-8 or abs(corners[0][1] - end0[1]) > 1e-8:
        return False, None, None, None

    v = []
    for s in segs:
        v_i = (float(s.p_end[0]) - float(s.p_start[0]), float(s.p_end[1]) - float(s.p_start[1]))
        v.append(v_i)

    def _dot(a, b): return a[0]*b[0] + a[1]*b[1]
    def _norm(a): return (a[0]**2 + a[1]**2) ** 0.5

    for i in range(4):
        j = (i + 1) % 4
        if abs(_dot(v[i], v[j])) > tol_orth * (_norm(v[i]) * _norm(v[j]) + 1e-16):
            return False, None, None, None

    horiz_idx, vert_idx = [], []
    for i in range(4):
        if abs(v[i][1]) <= abs(v[i][0]):
            horiz_idx.append(i)
        else:
            vert_idx.append(i)
    if len(horiz_idx) != 2 or len(vert_idx) != 2:
        return False, None, None, None

    def _len(idx): return _norm(v[idx])
    if abs(_len(horiz_idx[0]) - _len(horiz_idx[1])) > tol_len:
        return False, None, None, None
    if abs(_len(vert_idx[0]) - _len(vert_idx[1])) > tol_len:
        return False, None, None, None

    return True, corners, horiz_idx, vert_idx

# -----------------------------------------------------------------------------
# Build curves/loops and collect transfinite
# -----------------------------------------------------------------------------

def _build_curves_for_groups(
        gmsh_model, segment_groups: Sequence["SegmentGroup"], tolerance: float
) -> Tuple[List[int], List[int], List[Tuple[int, int]], Dict[int, Dict[str, object]], Dict[Tuple[float,float], int]]:
    """Create OCC curves for all groups.

    Returns
    -------
    closed_loops : list[int]
    internal_open : list[int]
    pending_transfinite : list[(curve_tag, n_subdiv)]
    loop_meta : dict loop_tag -> { group, is_rect, corners_pts, nx, ny }
    point_cache : dict key -> pointTag
    """
    occ = gmsh_model.occ
    pcache: Dict[Tuple[float, float], int] = {}
    closed_loops: List[int] = []
    internal_open: List[int] = []
    pending_transfinite: List[Tuple[int, int]] = []
    loop_meta: Dict[int, Dict[str, object]] = {}

    for group in segment_groups:
        segs = list(group.segments)
        curve_tags: List[int] = []
        for s in segs:
            ct, nsub = _add_curve_for_segment(occ, pcache, s)
            curve_tags.append(int(ct))
            if nsub is not None and int(nsub) >= 1:
                pending_transfinite.append((int(ct), int(nsub)))

        if group.is_closed(tolerance=tolerance):
            loop = int(occ.addCurveLoop(curve_tags))
            closed_loops.append(loop)

            is_rect, corner_coords, h_idx, v_idx = _is_rectangular_group(group)
            meta = {"group": group, "is_rect": is_rect, "corners_pts": None, "nx": None, "ny": None}
            if is_rect and corner_coords is not None:
                corner_tags: List[int] = []
                for c in corner_coords:
                    corner_tags.append(_add_point_cached(occ, pcache, c))
                def _n_of(i: int) -> Optional[int]:
                    n = getattr(segs[i], "n_subdivisions", None)
                    return int(n) if (n is not None and int(n) >= 1) else None
                nx_vals = [_n_of(h_idx[0]), _n_of(h_idx[1])]
                ny_vals = [_n_of(v_idx[0]), _n_of(v_idx[1])]
                nx = max([n for n in nx_vals if n is not None], default=None)
                ny = max([n for n in ny_vals if n is not None], default=None)
                meta.update({"corners_pts": list(map(int, corner_tags)), "nx": nx, "ny": ny})
            loop_meta[loop] = meta

        else:
            internal_open.extend(list(map(int, curve_tags)))

    return closed_loops, internal_open, pending_transfinite, loop_meta, pcache

# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------

def mesh_2d_gmsh(
        segment_groups: Sequence["SegmentGroup"],
        *,
        second_order: bool = False,
        mesh_type: int = 0,
        tolerance: float = 1e-6,
        name: str = "SegmentGroupMesh",
        open_fltk: bool = False,
) -> "Model":
    """Mesh SegmentGroups using native Gmsh geometry and return a Model.

    Notes
    -----
    - Transfinite for curves is applied AFTER occ.synchronize().
    - Perfect quad grid for rectangular outer boundary without holes.
    """
    import gmsh

    if not segment_groups:
        raise ValueError("mesh_2d_gmsh requires at least one SegmentGroup.")

    gmsh.initialize()
    try:
        gmsh.model.add(str(name) if name is not None else "SegmentGroupMesh")
        gmsh.option.setNumber("General.Terminal", 0)

        # Build geometry
        closed_loops, internal_curves, pending_transfinite, loop_meta, _pcache = _build_curves_for_groups(
            gmsh.model, segment_groups, tolerance
        )
        if not closed_loops:
            raise ValueError("At least one closed SegmentGroup is required to create a meshable surface.")

        outer_loop, *hole_loops = closed_loops
        surface_tag = int(gmsh.model.occ.addPlaneSurface([int(outer_loop)] + list(map(int, hole_loops))))

        # OCC → Model DB
        gmsh.model.occ.synchronize()

        # Transfinite curves now that curves exist in gmsh.model
        for curve_tag, nsub in pending_transfinite:
            gmsh.model.mesh.setTransfiniteCurve(int(curve_tag), int(nsub) + 1)

        # Perfect quad grid if rectangle and no holes
        meta_outer = loop_meta.get(int(outer_loop), {})
        is_rect = bool(meta_outer.get("is_rect", False))
        corners_pts = meta_outer.get("corners_pts")
        nx = meta_outer.get("nx")
        ny = meta_outer.get("ny")
        if is_rect and corners_pts and not hole_loops:
            # IMPORTANT: pass corners via keyword, otherwise 2nd positional arg is 'arrangement' (str)
            gmsh.model.mesh.setTransfiniteSurface(int(surface_tag), cornerTags=list(map(int, corners_pts)))
            gmsh.model.mesh.setRecombine(2, int(surface_tag))
        else:
            if mesh_type in {1, 2}:
                gmsh.option.setNumber("Mesh.RecombineAll", 1)

        # Split with internal open curves, if any
        if internal_curves:
            gmsh.model.occ.fragment([(2, int(surface_tag))], [(1, int(t)) for t in internal_curves])
            gmsh.model.occ.synchronize()

        # Meshing options
        gmsh.option.setNumber("Mesh.ElementOrder", 2 if second_order else 1)
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay

        # Generate mesh
        gmsh.model.mesh.generate(2)

        # Optional viewer
        if open_fltk:
            try:
                gmsh.fltk.run()
            except Exception:
                # headless / no FLTK available — ignore
                pass

        # Collect mesh → Model
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        coords = node_coords.reshape((-1, 3))

        from .model import Model as _Model  # late import avoids circular import
        mesh_model = _Model(str(name) if name is not None else "SegmentGroupMesh")

        tag_to_node: Dict[int, Node] = {}
        for tag, xyz in zip(node_tags, coords):
            x, y = float(xyz[0]), float(xyz[1])
            n = mesh_model.add_node(Node(None, x, y, 0.0))
            tag_to_node[int(tag)] = n

        # Elements → C2D*
        elt_types, _, elt_node_tags = gmsh.model.mesh.getElements(dim=2)
        gmsh_nodes_per_type = {2: 3, 3: 4, 9: 6, 10: 9}  # triangle3, quad4, tri6, quad9
        c2d_lookup = {3: "C2D3", 4: "C2D4", 6: "C2D6", 8: "C2D8"}

        for etype, flat_nodes in zip(elt_types, elt_node_tags):
            count = gmsh_nodes_per_type.get(int(etype))
            if count is None:
                continue
            arr = np.array(flat_nodes, dtype=int).reshape(-1, count)
            for row in arr:
                # map to Node objects; drop quad center node for 9-node to get C2D8
                node_objs = [tag_to_node[int(t)] for t in row]
                if count == 9:
                    node_objs = node_objs[:8]
                elem_type = c2d_lookup[len(node_objs)]
                ids = [int(n.node_id) for n in node_objs if n.node_id is not None]
                cls = Elements.element_class_for(elem_type)
                mesh_model.add_element(cls(None, ids))

        # NodeSets per group & per segment (polyline-based contains())
        group_sets: Dict[int, NodeSet] = {}
        segment_sets: Dict[int, NodeSet] = {}
        for group in segment_groups:
            gset = NodeSet(group.name, [])
            mesh_model.add_nodeset(gset)
            group_sets[id(group)] = gset
            for seg in group.segments:
                sset = NodeSet(seg.name, [])
                mesh_model.add_nodeset(sset)
                segment_sets[id(seg)] = sset

        for node in mesh_model.nodes._items:
            if node is None:
                continue
            p = (float(node.x), float(node.y))
            for group in segment_groups:
                gset = group_sets[id(group)]
                for seg in group.segments:
                    if seg.contains(np.asarray(p), tolerance=tolerance):
                        segment_sets[id(seg)].add(node)
                        gset.add(node)

        return mesh_model

    finally:
        if gmsh.isInitialized():
            gmsh.finalize()
