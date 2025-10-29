# file: run_buckle_with_circular_bspline_hole.py
from __future__ import annotations

import math
import random
import argparse
import tqdm
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import set_start_method
import sys
import yaml
from typing import List, Tuple, Dict, Iterable


if __package__ in (None, ""):
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    __package__ = "fembridge"

# ---- fembridge imports -------------------------------------------------------
from fembridge.core import Model, Runner
from fembridge.geometry import StraightSegment, SegmentGroup, BSplineSegment
from fembridge.materials import Material, ElasticityIsotropic, Density
from fembridge.sections.solid_section import SolidSection
from fembridge.sets.elementset import ElementSet
from fembridge.sets.nodeset import NodeSet
from fembridge.supports import Support, SupportCollector
from fembridge.loads import CLoad, LoadCollector
from fembridge.steps import LinearBucklingStep, LinearStaticStep  # ⟵ NEU


# =============================================================================
# Helpers for NodeSets
# =============================================================================
def _get_or_create_named_nodeset(model: Model, name: str, ns_or_iter) -> NodeSet:
    try:
        return model.node_sets[name]
    except KeyError:
        base = ns_or_iter if isinstance(ns_or_iter, NodeSet) else NodeSet.internal(ns_or_iter)
        created = NodeSet.create(name, base)
        model.add_nodeset(created)
        return created


# =============================================================================
# Tangential edge load builder (corners get half the weight)
# =============================================================================
def make_tangential_edge_loads(
        *,
        model: Model,
        width: float,
        height: float,
        ns_bottom: NodeSet,
        ns_right: NodeSet,
        ns_top: NodeSet,
        ns_left: NodeSet,
        name: str = "EDGE_TANGENTIAL",
) -> LoadCollector:
    def _count(ns: NodeSet) -> int:
        c = 0
        for _ in ns:
            c += 1
        return max(1, c)

    loads = LoadCollector(name)

    ns_bl = _get_or_create_named_nodeset(model, "BOTTOM_LEFT", (ns_bottom & ns_left))
    ns_br = _get_or_create_named_nodeset(model, "BOTTOM_RIGHT", (ns_bottom & ns_right))
    ns_tr = _get_or_create_named_nodeset(model, "TOP_RIGHT",   (ns_top    & ns_right))
    ns_tl = _get_or_create_named_nodeset(model, "TOP_LEFT",    (ns_top    & ns_left))

    def _edge(name_prefix: str, ns_edge: NodeSet, L: float, vec: tuple[float, float, float], cornerA: NodeSet, cornerB: NodeSet) -> None:
        N = _count(ns_edge)
        if N == 1:
            only = _get_or_create_named_nodeset(model, f"{name_prefix}_ALL", ns_edge)
            loads.add(CLoad(only, (vec[0]*L, vec[1]*L, vec[2]*L)))
            return

        w = float(L) / float(N - 1)
        w_end = 0.5 * w

        ns_corners_raw = (cornerA | cornerB)
        ns_corners = _get_or_create_named_nodeset(model, f"{name_prefix}_CORNERS", ns_corners_raw)

        ns_interior_raw = ns_edge - ns_corners
        ns_interior = _get_or_create_named_nodeset(model, f"{name_prefix}_INTERIOR", ns_interior_raw)

        if len(ns_interior) > 0:
            loads.add(CLoad(ns_interior, (vec[0]*w, vec[1]*w, vec[2]*w)))
        if len(ns_corners) > 0:
            loads.add(CLoad(ns_corners, (vec[0]*w_end, vec[1]*w_end, vec[2]*w_end)))

    _edge("BOTTOM", ns_bottom, float(width),  (+1.0, 0.0, 0.0), ns_bl, ns_br)
    _edge("RIGHT",  ns_right,  float(height), (0.0, -1.0, 0.0), ns_br, ns_tr)
    _edge("TOP",    ns_top,    float(width),  (-1.0, 0.0, 0.0), ns_tl, ns_tr)
    _edge("LEFT",   ns_left,   float(height), (0.0, +1.0, 0.0), ns_bl, ns_tl)

    return loads


# =============================================================================
# Hole generator under edge-clearance constraints
# =============================================================================
def _sample_hole_params_within_plate(
        *,
        width: float,
        height: float,
        radius_min: float,
        radius_max: float,
        jitter_xy: float,
        edge_clearance_min: float,
        n_ctrl_min: int,
        n_ctrl_max: int,
) -> tuple[float, float, float, int]:
    half_w = width / 2.0
    half_h = height / 2.0

    if radius_max <= 0:
        raise ValueError("radius_max must be > 0")
    if radius_min <= 0 or radius_min > radius_max:
        raise ValueError("radius_min must be in (0, radius_max]")
    if edge_clearance_min < 0:
        raise ValueError("edge_clearance_min must be >= 0")
    if n_ctrl_min < 3 or n_ctrl_max < n_ctrl_min:
        raise ValueError("n_ctrl_min >= 3 and n_ctrl_max >= n_ctrl_min")

    safety = radius_max + jitter_xy + edge_clearance_min
    if safety > min(half_w, half_h):
        raise ValueError(
            f"Infeasible hole constraints: radius_max({radius_max}) + jitter({jitter_xy}) + "
            f"edge_clearance_min({edge_clearance_min}) exceeds half of min(width, height)."
        )

    r = random.uniform(radius_min, radius_max)
    n_ctrl = random.randint(n_ctrl_min, n_ctrl_max)

    xmin = -half_w + (r + jitter_xy + edge_clearance_min)
    xmax = +half_w - (r + jitter_xy + edge_clearance_min)
    ymin = -half_h + (r + jitter_xy + edge_clearance_min)
    ymax = +half_h - (r + jitter_xy + edge_clearance_min)

    cx = random.uniform(xmin, xmax)
    cy = random.uniform(ymin, ymax)

    return cx, cy, r, n_ctrl


def _make_bspline_hole(
        *,
        cx: float,
        cy: float,
        radius: float,
        jitter_xy: float,
        n_ctrl: int,
        name: str = "HOLE",
) -> tuple[BSplineSegment, list[tuple[float, float]]]:
    angle_step = 2.0 * math.pi / n_ctrl
    ctrl: list[tuple[float, float]] = []
    for i in range(n_ctrl):
        a = i * angle_step
        x = cx + radius * math.cos(a) + random.uniform(-jitter_xy, +jitter_xy)
        y = cy + radius * math.sin(a) + random.uniform(-jitter_xy, +jitter_xy)
        ctrl.append((x, y))

    bspline = BSplineSegment(
        control_points=ctrl,
        closed=True,
        n_subdivisions=30,
        name=name,
    )
    return bspline, ctrl


# =============================================================================
# Curvature helpers (on FE hole nodes)
# =============================================================================
def _order_nodes_around_center(nodes_xy: List[Tuple[float, float]], center: Tuple[float, float]) -> List[int]:
    cx, cy = center
    angles = [math.atan2(y - cy, x - cx) for (x, y) in nodes_xy]
    order = sorted(range(len(nodes_xy)), key=lambda i: angles[i])
    return order

def _triangle_curvature(p0, p1, p2) -> float:
    # curvature k = 4*A / (|a||b||c|), where A is triangle area, signless here
    ax, ay = p0; bx, by = p1; cx, cy = p2
    a = math.hypot(bx - ax, by - ay)
    b = math.hypot(cx - bx, cy - by)
    c = math.hypot(ax - cx, ay - cy)
    if a == 0 or b == 0 or c == 0:
        return 0.0
    area2 = abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))  # 2*Area
    return (2.0 * area2) / (a * b * c)

def compute_curvature_on_hole_nodes(nodes_xy: List[Tuple[float, float]], center_xy: Tuple[float, float]) -> List[float]:
    if len(nodes_xy) < 3:
        return [0.0] * len(nodes_xy)
    order = _order_nodes_around_center(nodes_xy, center_xy)
    N = len(nodes_xy)
    curv = [0.0] * N
    for k in range(N):
        i = order[k]
        ip = order[(k - 1) % N]
        inx = order[(k + 1) % N]
        curv[i] = _triangle_curvature(nodes_xy[ip], nodes_xy[i], nodes_xy[inx])
    return curv


# =============================================================================
# Model builder (no-hole uses 2D C2D4 mesh; with-hole as before, then extrude)
# =============================================================================
def build_model(
        width: int,
        height: int,
        no_holes: bool,
        elem_size: float = 2.0,
        *,
        edge_clearance_min: float = 4.0,
        radius_min: float = 4.0,
        radius_max: float = 12.0,
        jitter_xy: float = 2.0,
        n_ctrl_min: int = 3,
        n_ctrl_max: int = 5,
) -> tuple[Model, SupportCollector, LoadCollector, dict, dict]:
    half_w = width / 2.0
    half_h = height / 2.0

    nx = max(1, int(round(width  / elem_size)))
    ny = max(1, int(round(height / elem_size)))

    seg_bottom = StraightSegment((-half_w, -half_h), (+half_w, -half_h), n_subdivisions=nx, name="BOTTOM")
    seg_right  = StraightSegment((+half_w, -half_h), (+half_w, +half_h), n_subdivisions=ny, name="RIGHT")
    seg_top    = StraightSegment((+half_w, +half_h), (-half_w, +half_h), n_subdivisions=nx, name="TOP")
    seg_left   = StraightSegment((-half_w, +half_h), (-half_w, -half_h), n_subdivisions=ny, name="LEFT")
    outer = SegmentGroup([seg_bottom, seg_right, seg_top, seg_left], name="OUTER")

    boundaries = [outer]
    hole_meta = {"present": False}

    if not no_holes:
        cx, cy, radius, n_ctrl = _sample_hole_params_within_plate(
            width=width, height=height,
            radius_min=radius_min, radius_max=radius_max,
            jitter_xy=jitter_xy, edge_clearance_min=edge_clearance_min,
            n_ctrl_min=n_ctrl_min, n_ctrl_max=n_ctrl_max,
        )
        bspline, ctrl_pts = _make_bspline_hole(
            cx=cx, cy=cy, radius=radius, jitter_xy=jitter_xy, n_ctrl=n_ctrl, name="HOLE"
        )
        hole = SegmentGroup([bspline], name="HOLE")
        boundaries = [outer, hole]
        hole_meta = {
            "present": True,
            "center_xy": [float(cx), float(cy)],
            "radius": float(radius),
            "jitter_xy": float(jitter_xy),
            "n_ctrl": int(n_ctrl),
            "control_points_xy": [[float(x), float(y)] for (x, y) in ctrl_pts],
            "edge_clearance_min": float(edge_clearance_min),
            "note": "Uniform-angle ctrl points with ±jitter_xy; center sampled to satisfy clearance.",
        }

    plate2d = Model.mesh_2d(
        boundaries,
        mesh_type=2,
        name=f"Plate{width}x{height}{'_noHole' if no_holes else '_withHole'}_2D",
        open_fltk=False
    )

    solid = plate2d.extruded(n=5, spacing=0.2)
    solid.name = f"SolidPlate_{width}x{height}_t1p0_n5{'_noHole' if no_holes else '_bspline'}"

    all_elems = [e for e in solid.elements._items if e is not None]
    elem_set = ElementSet("EALL", all_elems)
    solid.add_elementset(elem_set)

    steel = solid.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e3, 0.30))
        .set_density(Density(7850.0e-12))
    )
    solid.add_section(SolidSection(material=steel, elset=elem_set))

    ns_bottom = solid.node_sets["BOTTOM"]
    ns_right  = solid.node_sets["RIGHT"]
    ns_top    = solid.node_sets["TOP"]
    ns_left   = solid.node_sets["LEFT"]

    ns_bl = _get_or_create_named_nodeset(solid, "BOTTOM_LEFT",  (ns_bottom & ns_left))
    ns_br = _get_or_create_named_nodeset(solid, "BOTTOM_RIGHT", (ns_bottom & ns_right))
    ns_tr = _get_or_create_named_nodeset(solid, "TOP_RIGHT",    (ns_top    & ns_right))
    ns_tl = _get_or_create_named_nodeset(solid, "TOP_LEFT",     (ns_top    & ns_left))

    supports = SupportCollector("SUPPORTS")
    supports.add(Support(ns_bl, (0.0, 0.0, None, None, None, None)))
    supports.add(Support(ns_br, (None, 0.0, None, None, None, None)))
    for edge in (ns_bottom, ns_right, ns_top, ns_left):
        supports.add(Support(edge, (None, None, 0.0, None, None, None)))
    solid.add_supportcollector(supports)

    loads = make_tangential_edge_loads(
        model=solid,
        width=width,
        height=height,
        ns_bottom=ns_bottom,
        ns_right=ns_right,
        ns_top=ns_top,
        ns_left=ns_left,
        name="EDGE_TANGENTIAL",
    )
    solid.add_loadcollector(loads)

    meta = {
        "width": int(width),
        "height": int(height),
        "plate_centered": True,
        "mesh_type": 2,
        "mesh_target_edge_size": float(elem_size),
        "nx": int(nx),
        "ny": int(ny),
        "extruded_layers": 5,
        "layer_spacing": 0.2,
        "hole_params": {
            "used": not no_holes,
            "edge_clearance_min": float(edge_clearance_min),
            "radius_min": float(radius_min),
            "radius_max": float(radius_max),
            "jitter_xy": float(jitter_xy),
            "n_ctrl_min": int(n_ctrl_min),
            "n_ctrl_max": int(n_ctrl_max),
        },
    }
    return solid, supports, loads, hole_meta, meta


# =============================================================================
# Solution helpers
# =============================================================================
def first_buckling_value(solution) -> float:
    return float(solution.steps[1 if len(solution.steps) > 1 else 0].fields["BUCKLING_FACTORS"][0].item())

def collect_hole_nodes_at_midplane(model: Model) -> List:
    try:
        ns_hole = model.node_sets["HOLE"]
    except KeyError:
        return []
    out_nodes = []
    seen = set()
    for n in ns_hole:
        if abs(float(n.z) - 0.0) < 1e-12:
            key = (int(n.node_id), round(float(n.x), 12), round(float(n.y), 12))
            if key not in seen:
                seen.add(key)
                out_nodes.append(n)
    return out_nodes

def _mises_from_voigt6(s: Iterable[float]) -> float:
    # s = [sx, sy, sz, sxy, syz, szx]  (VTK order!)
    sx, sy, sz, sxy, syz, szx = map(float, s)
    sxx, syy, szz = sx, sy, sz
    txy, tyz, tzx = sxy, syz, szx
    term = ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2 + 6.0 * (txy ** 2 + tyz ** 2 + tzx ** 2))
    return math.sqrt(0.5 * term)

def extract_node_fields(solution, step_index: int, node_ids: List[int], want_fields: List[str]) -> Dict[str, List]:
    """Extract selected nodal fields by name (case-insensitive contains)."""
    out: Dict[str, List] = {}
    fields = solution.steps[step_index].fields
    # map candidate names to actual
    for wf in want_fields:
        actual = next((k for k in fields.keys() if wf.lower() in k.lower()), None)
        if actual is None:
            continue
        arr = fields[actual]
        # assume arr is (N, m) or (N,)
        if hasattr(arr, "__len__"):
            vals = []
            for nid in node_ids:
                vals.append(arr[nid] if nid < len(arr) else None)
            out[actual] = vals
    # add MISES if 6-comp stress found
    stress_key = next((k for k in out.keys() if "stress" in k.lower()), None)
    if stress_key is not None:
        S = out[stress_key]
        try:
            mises = []
            for row in S:
                if row is None:
                    mises.append(None)
                else:
                    # ensure order is VTK: XX,YY,ZZ,XY,YZ,XZ ; FEMaster array may already be reordered upstream
                    vv = list(map(float, row))
                    if len(vv) == 6:
                        mises.append(_mises_from_voigt6(vv))
                    else:
                        mises.append(None)
            out["MISES"] = mises
        except Exception:
            pass
    return out


# =============================================================================
# One run → YAML path
# =============================================================================
def run_single(
        out_dir: Path,
        engine_path: Path,
        width: int,
        height: int,
        no_holes: bool,
        elem_size: float,
        *,
        edge_clearance_min: float,
        radius_min: float,
        radius_max: float,
        jitter_xy: float,
        n_ctrl_min: int,
        n_ctrl_max: int,
        compute_stresses: bool,
        compute_curvature: bool,
) -> Path:
    random.seed()

    model, supports, loads, hole_meta, plate_meta = build_model(
        width=width,
        height=height,
        no_holes=no_holes,
        elem_size=elem_size,
        edge_clearance_min=edge_clearance_min,
        radius_min=radius_min,
        radius_max=radius_max,
        jitter_xy=jitter_xy,
        n_ctrl_min=n_ctrl_min,
        n_ctrl_max=n_ctrl_max,
    )

    # Steps:
    #  - optional LinearStaticStep (for stresses)
    #  - LinearBucklingStep (as before)
    steps = []
    if compute_stresses:
        steps.append(LinearStaticStep(load_collectors=[loads], support_collectors=[supports]))
    steps.append(LinearBucklingStep(sigma=1, num_eigenvalues=1, load_collectors=[loads], support_collectors=[supports]))
    for st in steps:
        model.add_step(st)

    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
    runner.set_option(Runner.Option.NO_TEMP_FILES, False)
    solution = runner.run()

    lam1 = first_buckling_value(solution)

    # hole nodes @ midplane with IDs and XY
    hole_nodes = collect_hole_nodes_at_midplane(model)
    hole_nodes_xy = [(float(n.x), float(n.y)) for n in hole_nodes]
    hole_node_ids = [int(n.node_id) for n in hole_nodes]

    curvature_vals: List[float] = []
    if compute_curvature and hole_meta.get("present", False):
        cx, cy = hole_meta["center_xy"]
        curvature_vals = compute_curvature_on_hole_nodes(hole_nodes_xy, (cx, cy))

    node_outputs: Dict[str, List] = {}
    if compute_stresses and len(solution.steps) >= 1:
        static_step_index = 0  # we appended static first
        node_outputs = extract_node_fields(
            solution,
            static_step_index,
            hole_node_ids,
            want_fields=["STRESS", "DISPLACEMENT"]  # add more if you like
        )

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    out_dir.mkdir(parents=True, exist_ok=True)
    hole_tag = "nohole" if no_holes else "bspline"
    out_path = out_dir / f"buckling_{hole_tag}_w{width}_h{height}_{timestamp}.yaml"

    # pack per-node table
    per_node = []
    for i, n in enumerate(hole_nodes):
        entry = {
            "node_id": int(n.node_id),
            "x": float(n.x),
            "y": float(n.y),
            "z": float(n.z),
        }
        if compute_curvature and curvature_vals:
            entry["curvature"] = float(curvature_vals[i])
        # attach any extracted node fields
        for k, arr in node_outputs.items():
            val = arr[i]
            if hasattr(val, "__len__") and not isinstance(val, (str, bytes)):
                entry[k] = [float(v) for v in val]  # cast to float list
            elif val is None:
                entry[k] = None
            else:
                entry[k] = float(val)
        per_node.append(entry)

    payload = {
        "meta": {
            "timestamp": timestamp,
            "plate": {
                "min_xy": [-plate_meta["width"]/2.0, -plate_meta["height"]/2.0],
                "max_xy": [ plate_meta["width"]/2.0,  plate_meta["height"]/2.0],
                "width": int(plate_meta["width"]),
                "height": int(plate_meta["height"]),
                "extruded_layers": plate_meta.get("extruded_layers", 5),
                "layer_spacing": plate_meta.get("layer_spacing", 0.2),
            },
            "material": {"E_MPa": 210e3, "nu": 0.30, "rho_t_per_mm3": 7850.0e-12},
            "loads": "outer edges tangential (+x bottom, -y right, -x top, +y left), sum per edge = edge length; corner nodes get half",
            "supports": "BL: u1=u2=0; BR: u2=0; all outer edges: w=0",
            "mesh": {
                "type": plate_meta["mesh_type"],
                "target_edge_size": plate_meta["mesh_target_edge_size"],
                "nx": plate_meta["nx"],
                "ny": plate_meta["ny"],
            },
            "options": {
                "compute_stresses": bool(compute_stresses),
                "compute_curvature": bool(compute_curvature),
            }
        },
        "hole": hole_meta,
        "first_buckling_factor": float(lam1),
        "hole_nodes": per_node,  # ⟵ Spannungen/Krümmung pro HOLE-Knoten (Midplane)
    }

    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

    return out_path


def _run_once(
        seed: int,
        out_dir: Path,
        engine_path: Path,
        width: int,
        height: int,
        no_holes: bool,
        elem_size: float,
        *,
        edge_clearance_min: float,
        radius_min: float,
        radius_max: float,
        jitter_xy: float,
        n_ctrl_min: int,
        n_ctrl_max: int,
        compute_stresses: bool,
        compute_curvature: bool,
):
    random.seed(seed)
    try:
        p = run_single(
            out_dir, engine_path, width=width, height=height, no_holes=no_holes, elem_size=elem_size,
            edge_clearance_min=edge_clearance_min, radius_min=radius_min, radius_max=radius_max,
            jitter_xy=jitter_xy, n_ctrl_min=n_ctrl_min, n_ctrl_max=n_ctrl_max,
            compute_stresses=compute_stresses, compute_curvature=compute_curvature
        )
        return True, seed, (width, height), str(p)
    except Exception as e:
        return False, seed, (width, height), str(e)


# =============================================================================
# CLI helpers
# =============================================================================
def parse_seq(spec: str, default_step: int = 20) -> List[int]:
    parts = spec.split(":")
    if len(parts) not in (2, 3):
        raise ValueError("range must be 'start:end[:step]'")
    start = int(parts[0]); end = int(parts[1])
    step = int(parts[2]) if len(parts) == 3 else default_step
    if step <= 0:
        raise ValueError("step must be > 0")
    if end < start:
        raise ValueError("end must be >= start")
    vals = list(range(start, end + 1, step))
    if vals[-1] != end:
        vals.append(end)
    return vals


# =============================================================================
# CLI
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Generate buckling samples (2D mesh → extrude). Optional: linear static stresses on HOLE nodes and curvature."
    )
    parser.add_argument("--n", type=int, required=True,
                        help="Number of samples to generate (per width/height combo if seq is used).")
    parser.add_argument("--out", type=Path, required=True, help="Output directory for YAML files.")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers.")
    parser.add_argument("--seed", type=int, default=None, help="Base RNG seed (optional).")
    parser.add_argument("--engine", type=Path, default=None, help="Path to FEMaster engine (overrides default).")

    parser.add_argument("--no_holes", action="store_true", help="Use rectangle-only 2D mesh (C2D4), then extrude.")
    parser.add_argument("--width", type=int, default=40, help="Plate width (ignored if --width_seq is given).")
    parser.add_argument("--height", type=int, default=40, help="Plate height (ignored if --height_seq is given).")
    parser.add_argument("--width_seq", type=str, default=None,
                        help="Generate for multiple widths: 'start:end[:step]'. --n is per (w,h).")
    parser.add_argument("--height_seq", type=str, default=None,
                        help="Generate for multiple heights: 'start:end[:step]'. --n is per (w,h).")
    parser.add_argument("--elem_size", type=float, default=2.0, help="Target in-plane element size for meshing.")

    # Hole parameter controls
    parser.add_argument("--edge_clearance_min", type=float, default=4.0,
                        help="Minimal distance from hole (radius+jitter) to each plate edge.")
    parser.add_argument("--radius_min", type=float, default=4.0, help="Minimum hole radius.")
    parser.add_argument("--radius_max", type=float, default=12.0, help="Maximum hole radius.")
    parser.add_argument("--jitter_xy", type=float, default=2.0,
                        help="±Jitter applied equally in x and y to each control point.")
    parser.add_argument("--n_ctrl_min", type=int, default=3, help="Minimum number of B-spline control points (>=3).")
    parser.add_argument("--n_ctrl_max", type=int, default=5, help="Maximum number of B-spline control points.")

    # Optional analyses on HOLE nodes
    parser.add_argument("--compute_stresses", action="store_true",
                        help="Add a LinearStatic step (using same loads/supports) and extract nodal fields on HOLE nodes.")
    parser.add_argument("--compute_curvature", action="store_true",
                        help="Compute curvature on HOLE nodes (midplane) from their polygonal geometry.")

    args = parser.parse_args()

    default_engine = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    engine_path = args.engine if args.engine is not None else default_engine

    try:
        set_start_method("spawn")
    except RuntimeError:
        pass

    widths = parse_seq(args.width_seq) if args.width_seq else [int(args.width)]
    heights = parse_seq(args.height_seq) if args.height_seq else [int(args.height)]

    base_rng = random.Random(args.seed) if args.seed is not None else random.Random()
    total_runs = args.n * (len(widths) * len(heights))
    seeds: List[int] = [base_rng.getrandbits(64) for _ in range(total_runs)]

    ok = fail = 0
    args.out.mkdir(parents=True, exist_ok=True)

    jobs: List[Tuple[int, int]] = []
    for w in widths:
        for h in heights:
            for _ in range(args.n):
                jobs.append((w, h))

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for i, (w, h) in enumerate(jobs):
            futs.append(ex.submit(
                _run_once, seeds[i], args.out, engine_path, w, h, args.no_holes, args.elem_size,
                edge_clearance_min=args.edge_clearance_min,
                radius_min=args.radius_min,
                radius_max=args.radius_max,
                jitter_xy=args.jitter_xy,
                n_ctrl_min=args.n_ctrl_min,
                n_ctrl_max=args.n_ctrl_max,
                compute_stresses=args.compute_stresses,
                compute_curvature=args.compute_curvature,
            ))

        with tqdm.tqdm(total=len(futs), desc="Buckling runs", ncols=100, file=sys.stdout) as bar:
            for fut in as_completed(futs):
                success, seed, wh, msg = fut.result()
                if success:
                    ok += 1
                else:
                    fail += 1
                    tqdm.tqdm.write(f"[seed={seed} w={wh[0]} h={wh[1]}] ERROR: {msg}", file=sys.stdout)
                bar.update(1)
                bar.set_postfix(ok=ok, fail=fail)

    print(f"\nDone. ok={ok}, failed={fail}. Output dir: {args.out.resolve()}")


if __name__ == "__main__":
    main()
