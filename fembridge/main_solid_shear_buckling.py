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
from typing import List, Tuple


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
from fembridge.steps import LinearBucklingStep


# =============================================================================
# Abaqus-like random B-spline hole (for the with-hole branch)
# =============================================================================
def make_abaqus_like_bspline(
        *,
        radius_min: float = 4.0,
        radius_max: float = 12.0,
        jitter_xy: float = 2.0,
        n_ctrl_min: int = 3,
        n_ctrl_max: int = 5,
        name: str = "HOLE",
) -> tuple[BSplineSegment, list[tuple[float, float]], float]:
    radius = random.uniform(radius_min, radius_max)
    n_ctrl = random.randint(n_ctrl_min, n_ctrl_max)

    angle_step = 2.0 * math.pi / n_ctrl
    ctrl: list[tuple[float, float]] = []
    for i in range(n_ctrl):
        a = i * angle_step
        x = radius * math.cos(a) + random.uniform(-jitter_xy, +jitter_xy)
        y = radius * math.sin(a) + random.uniform(-jitter_xy, +jitter_xy)
        ctrl.append((x, y))

    bspline = BSplineSegment(
        control_points=ctrl,
        closed=True,
        n_subdivisions=30,
        name=name,
    )
    return bspline, ctrl, radius


# =============================================================================
# Helpers for NodeSets
# =============================================================================
def _get_or_create_named_nodeset(model: Model, name: str, ns_or_iter) -> NodeSet:
    """
    Fetch a NodeSet by name if it exists; otherwise create it from `ns_or_iter`
    and add it to the model with the given name.
    """
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
    """
    Creates tangential edge loads on the four outer edges with a trapezoidal
    node distribution so that:
      - interior nodes get w = L/(N-1)
      - end (corner) nodes get w/2
      - sum over a given edge equals its geometric length L

    Directions:
      bottom: +x, right: -y, top: -x, left: +y

    All intermediate NodeSets are **named** and added to the model to avoid
    'internal set' semantics (which would otherwise take only the first id).
    """

    def _count(ns: NodeSet) -> int:
        c = 0
        for _ in ns:
            c += 1
        return max(1, c)

    loads = LoadCollector(name)

    # Corner sets (re-use if already present, else create & add)
    ns_bl = _get_or_create_named_nodeset(model, "BOTTOM_LEFT", (ns_bottom & ns_left))
    ns_br = _get_or_create_named_nodeset(model, "BOTTOM_RIGHT", (ns_bottom & ns_right))
    ns_tr = _get_or_create_named_nodeset(model, "TOP_RIGHT",   (ns_top    & ns_right))
    ns_tl = _get_or_create_named_nodeset(model, "TOP_LEFT",    (ns_top    & ns_left))

    # For each edge: create named CORNERS and INTERIOR sets, add to model, then place loads
    def _edge(name_prefix: str, ns_edge: NodeSet, L: float, vec: tuple[float, float, float], cornerA: NodeSet, cornerB: NodeSet) -> None:
        N = _count(ns_edge)
        if N == 1:
            # Degenerate: put full length onto the single node
            only = _get_or_create_named_nodeset(model, f"{name_prefix}_ALL", ns_edge)
            loads.add(CLoad(only, (vec[0]*L, vec[1]*L, vec[2]*L)))
            return

        w = float(L) / float(N - 1)  # interior node weight
        w_end = 0.5 * w              # end node weight

        ns_corners_raw = (cornerA | cornerB)
        ns_corners = _get_or_create_named_nodeset(model, f"{name_prefix}_CORNERS", ns_corners_raw)

        ns_interior_raw = ns_edge - ns_corners
        ns_interior = _get_or_create_named_nodeset(model, f"{name_prefix}_INTERIOR", ns_interior_raw)

        if len(ns_interior) > 0:
            loads.add(CLoad(ns_interior, (vec[0]*w, vec[1]*w, vec[2]*w)))
        if len(ns_corners) > 0:
            loads.add(CLoad(ns_corners, (vec[0]*w_end, vec[1]*w_end, vec[2]*w_end)))

    # bottom: +x
    _edge("BOTTOM", ns_bottom, float(width),  (+1.0, 0.0, 0.0), ns_bl, ns_br)
    # right:  -y
    _edge("RIGHT",  ns_right,  float(height), (0.0, -1.0, 0.0), ns_br, ns_tr)
    # top:    -x
    _edge("TOP",    ns_top,    float(width),  (-1.0, 0.0, 0.0), ns_tl, ns_tr)
    # left:   +y
    _edge("LEFT",   ns_left,   float(height), (0.0, +1.0, 0.0), ns_bl, ns_tl)

    return loads


# =============================================================================
# Model builder (no-hole uses 2D C2D4 mesh; with-hole as before, then extrude)
# =============================================================================
def build_model(
        width: int,
        height: int,
        no_holes: bool,
        elem_size: float = 2.0,
) -> tuple[Model, SupportCollector, LoadCollector, dict, dict]:
    """
    - no_holes=True  → 2D-Rect (C2D4) ohne Innenkontur, dann extrudiert.
    - no_holes=False → Outer + HOLE-Loop, 2D (C2D4/C2D3 je nach Topologie), dann extrudiert.
    """

    half_w = width / 2.0
    half_h = height / 2.0

    # discrete subdivisions derived from target elem_size
    nx = max(1, int(round(width  / elem_size)))
    ny = max(1, int(round(height / elem_size)))

    # Outer boundary (CCW, start at bottom edge)
    seg_bottom = StraightSegment((-half_w, -half_h), (+half_w, -half_h), n_subdivisions=nx, name="BOTTOM")
    seg_right  = StraightSegment((+half_w, -half_h), (+half_w, +half_h), n_subdivisions=ny, name="RIGHT")
    seg_top    = StraightSegment((+half_w, +half_h), (-half_w, +half_h), n_subdivisions=nx, name="TOP")
    seg_left   = StraightSegment((-half_w, +half_h), (-half_w, -half_h), n_subdivisions=ny, name="LEFT")
    outer = SegmentGroup([seg_bottom, seg_right, seg_top, seg_left], name="OUTER")

    boundaries = [outer]
    hole_meta = {"present": False}

    if not no_holes:
        bspline, ctrl_pts, radius = make_abaqus_like_bspline(name="HOLE")
        hole = SegmentGroup([bspline], name="HOLE")
        boundaries = [outer, hole]
        hole_meta = {
            "present": True,
            "radius": float(radius),
            "control_points_xy": [[float(x), float(y)] for (x, y) in ctrl_pts],
            "note": "Abaqus-like: uniform angles with ±2 x/y jitter",
        }

    # --- 2D meshing (C2D4/C2D3 depending on topology), then extrude ---
    plate2d = Model.mesh_2d(
        boundaries,
        mesh_type=2,  # quads where possible (rectangles get transfinite quads)
        name=f"Plate{width}x{height}{'_noHole' if no_holes else '_withHole'}_2D",
        open_fltk=False
    )

    solid = plate2d.extruded(n=5, spacing=0.2)
    solid.name = f"SolidPlate_{width}x{height}_t1p0_n5{'_noHole' if no_holes else '_bspline'}"

    # --- ElementSet + Section/Material ---
    all_elems = [e for e in solid.elements._items if e is not None]
    elem_set = ElementSet("EALL", all_elems)
    solid.add_elementset(elem_set)

    steel = solid.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e3, 0.30))   # MPa (210 GPa)
        .set_density(Density(7850.0e-12))                  # tonne/mm^3 (if mm units)
    )
    solid.add_section(SolidSection(material=steel, elset=elem_set))

    # --- NodeSets from mesher (names propagated from segments) ---
    ns_bottom = solid.node_sets["BOTTOM"]
    ns_right  = solid.node_sets["RIGHT"]
    ns_top    = solid.node_sets["TOP"]
    ns_left   = solid.node_sets["LEFT"]

    # Corner sets (make them persistent with names; if they already exist, getter returns them)
    ns_bl = _get_or_create_named_nodeset(solid, "BOTTOM_LEFT",  (ns_bottom & ns_left))
    ns_br = _get_or_create_named_nodeset(solid, "BOTTOM_RIGHT", (ns_bottom & ns_right))
    ns_tr = _get_or_create_named_nodeset(solid, "TOP_RIGHT",    (ns_top    & ns_right))
    ns_tl = _get_or_create_named_nodeset(solid, "TOP_LEFT",     (ns_top    & ns_left))

    # --- Supports ---
    supports = SupportCollector("SUPPORTS")
    supports.add(Support(ns_bl, (0.0, 0.0, None, None, None, None)))   # BL: u1=0, u2=0
    supports.add(Support(ns_br, (None, 0.0, None, None, None, None)))  # BR:      u2=0
    for edge in (ns_bottom, ns_right, ns_top, ns_left):
        supports.add(Support(edge, (None, None, 0.0, None, None, None)))  # w=0
    solid.add_supportcollector(supports)

    # --- Edge loads (corners half, sum per edge equals edge length) ---
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
    }
    return solid, supports, loads, hole_meta, meta


# =============================================================================
# Solution helpers
# =============================================================================
def first_buckling_value(solution) -> float:
    return float(solution.steps[0].fields["BUCKLING_FACTORS"][0].item())


def collect_hole_nodes_xy_at_midplane(model: Model) -> list[tuple[float, float]]:
    # If there is no HOLE set → empty list
    try:
        ns_hole = model.node_sets["HOLE"]
    except KeyError:
        return []
    out_xy: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for n in ns_hole:
        if abs(float(n.z) - 0.0) < 1e-12:
            xy = (round(float(n.x), 12), round(float(n.y), 12))
            if xy not in seen:
                seen.add(xy)
                out_xy.append((float(n.x), float(n.y)))
    return out_xy


# =============================================================================
# One run → YAML path
# =============================================================================
def run_single(out_dir: Path, engine_path: Path, width: int, height: int, no_holes: bool, elem_size: float) -> Path:
    random.seed()

    model, supports, loads, hole_meta, plate_meta = build_model(
        width=width, height=height, no_holes=no_holes, elem_size=elem_size
    )
    model.add_step(
        LinearBucklingStep(
            sigma=1,
            num_eigenvalues=1,
            load_collectors=[loads],
            support_collectors=[supports],
        )
    )

    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
    runner.set_option(Runner.Option.NO_TEMP_FILES, False)
    solution = runner.run()

    lam1 = first_buckling_value(solution)
    hole_xy_mid = collect_hole_nodes_xy_at_midplane(model)

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    out_dir.mkdir(parents=True, exist_ok=True)
    hole_tag = "nohole" if no_holes else "bspline"
    out_path = out_dir / f"buckling_{hole_tag}_w{width}_h{height}_{timestamp}.yaml"

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
            }
        },
        "hole": hole_meta,
        "first_buckling_factor": float(lam1),
        "hole_nodes_xy": [[x, y] for (x, y) in hole_xy_mid],
    }
    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

    return out_path


def _run_once(seed: int, out_dir: Path, engine_path: Path, width: int, height: int, no_holes: bool, elem_size: float):
    random.seed(seed)
    try:
        p = run_single(out_dir, engine_path, width=width, height=height, no_holes=no_holes, elem_size=elem_size)
        return True, seed, (width, height), str(p)
    except Exception as e:
        return False, seed, (width, height), str(e)


# =============================================================================
# CLI helpers
# =============================================================================
def parse_seq(spec: str, default_step: int = 20) -> List[int]:
    """
    Parse "start:end[:step]" into an inclusive int list.
    Example: "40:200:20" -> [40,60,80,...,200]
             "40:200"    -> same with default step
    """
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
        description="Generate buckling samples (2D mesh → extrude; optional B-spline hole)."
    )
    parser.add_argument("--n", type=int, required=True,
                        help="Number of samples to generate (per width if width range is used).")
    parser.add_argument("--out", type=Path, required=True, help="Output directory for YAML files.")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers.")
    parser.add_argument("--seed", type=int, default=None, help="Base RNG seed (optional).")
    parser.add_argument("--engine", type=Path, default=None, help="Path to FEMaster engine (overrides default).")

    parser.add_argument("--no_holes", action="store_true", help="Use rectangle-only 2D mesh (C2D4), then extrude.")
    parser.add_argument("--width", type=int, default=40, help="Plate width (ignored if --width_seq is given).")
    parser.add_argument("--height", type=int, default=40, help="Plate height (fixed if only width varies).")
    parser.add_argument("--width_seq", type=str, default=None,
                        help="Generate for multiple widths: 'start:end[:step]'. --n is per width.")
    parser.add_argument("--elem_size", type=float, default=2.0, help="Target in-plane element size for meshing.")

    args = parser.parse_args()

    default_engine = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    engine_path = args.engine if args.engine is not None else default_engine

    try:
        set_start_method("spawn")
    except RuntimeError:
        pass

    # width list
    widths = parse_seq(args.width_seq) if args.width_seq else [int(args.width)]
    height = int(args.height)

    base_rng = random.Random(args.seed) if args.seed is not None else random.Random()
    total_runs = args.n * len(widths)
    seeds: List[int] = [base_rng.getrandbits(64) for _ in range(total_runs)]

    ok = fail = 0
    args.out.mkdir(parents=True, exist_ok=True)

    jobs: List[Tuple[int, int]] = []
    for w in widths:
        for _ in range(args.n):
            jobs.append((w, height))

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for i, (w, h) in enumerate(jobs):
            futs.append(ex.submit(_run_once, seeds[i], args.out, engine_path, w, h, args.no_holes, args.elem_size))

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
