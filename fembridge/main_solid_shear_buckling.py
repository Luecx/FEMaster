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

from fembridge.core import Model, Runner
from fembridge.geometry import StraightSegment, SegmentGroup, BSplineSegment
from fembridge.materials import Material, ElasticityIsotropic, Density
from fembridge.sections.solid_section import SolidSection
from fembridge.sets.elementset import ElementSet
from fembridge.supports import Support, SupportCollector
from fembridge.loads import CLoad, LoadCollector
from fembridge.steps import LinearBucklingStep


# -----------------------------------------------------------------------------
# Geometry helpers (Abaqus-style randomness) — only used if holes are enabled
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Model builder
# -----------------------------------------------------------------------------
def build_model(width: int, height: int, no_holes: bool) -> tuple[Model, SupportCollector, LoadCollector, dict]:
    """
    Erstellt das 3D-Modell (extrudierte Platte) mit/ohne Loch.
    WICHTIG: Nur die BREITE variiert; die HÖHE ist fix (über --height).
    """
    # Outer rectangle: centered at origin
    # Ziel-Kantenlänge ≈ 2.0 ⇒ Unterteilungen nx ~ width/2, ny ~ height/2
    nx = max(1, int(round(width / 2)))
    ny = max(1, int(round(height / 2)))
    half_w = width / 2.0
    half_h = height / 2.0

    seg_bottom = StraightSegment((-half_w, -half_h), (+half_w, -half_h), n_subdivisions=nx, name="BOTTOM")
    seg_right  = StraightSegment((+half_w, -half_h), (+half_w, +half_h), n_subdivisions=ny, name="RIGHT")
    seg_top    = StraightSegment((+half_w, +half_h), (-half_w, +half_h), n_subdivisions=nx, name="TOP")
    seg_left   = StraightSegment((-half_w, +half_h), (-half_w, -half_h), n_subdivisions=ny, name="LEFT")
    outer = SegmentGroup([seg_bottom, seg_right, seg_top, seg_left], name="OUTER")

    # Hole (optional)
    hole_meta = {"present": False}
    boundaries = [outer]
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

    # Mesh 2D (quads), then extrude to thickness 1.0 (5 layers × 0.2)
    plate2d = Model.mesh_2d(boundaries, mesh_type=2, name=f"Plate{width}x{height}{'_noHole' if no_holes else '_withHole'}_2D")
    solid = plate2d.extruded(n=5, spacing=0.2)
    solid.name = f"SolidPlate_{width}x{height}_t1p0_n5{'_noHole' if no_holes else '_bspline'}"

    # Elements & section
    all_elems = [e for e in solid.elements._items if e is not None]
    elem_set = ElementSet("EALL", all_elems)
    solid.add_elementset(elem_set)

    steel = solid.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e3, 0.30))   # MPa (210 GPa)
        .set_density(Density(7850.0e-12))                  # tonne/mm^3 (bei mm)
    )
    solid.add_section(SolidSection(material=steel, elset=elem_set))

    # Node sets (edges auto-named; persist corner sets)
    ns_bottom = solid.node_sets["BOTTOM"]
    ns_right  = solid.node_sets["RIGHT"]
    ns_top    = solid.node_sets["TOP"]
    ns_left   = solid.node_sets["LEFT"]

    ns_bl = ns_bottom & ns_left
    ns_br = ns_bottom & ns_right
    ns_bl.name = "BOTTOM_LEFT"
    ns_br.name = "BOTTOM_RIGHT"
    solid.add_nodeset(ns_bl)
    solid.add_nodeset(ns_br)

    # Supports
    supports = SupportCollector("SUPPORTS")
    supports.add(Support(ns_bl, (0.0, 0.0, None, None, None, None)))   # BL: u1=0, u2=0
    supports.add(Support(ns_br, (None, 0.0, None, None, None, None)))  # BR:      u2=0
    for edge in (ns_bottom, ns_right, ns_top, ns_left):
        supports.add(Support(edge, (None, None, 0.0, None, None, None)))  # w=0
    solid.add_supportcollector(supports)

    # Edge loads: Sum per edge equals its length (width for bottom/top, height for right/left).
    def count_nodes(ns) -> int:
        c = 0
        for _ in ns:
            c += 1
        return max(1, c)

    N_bot = count_nodes(ns_bottom)
    N_rig = count_nodes(ns_right)
    N_top = count_nodes(ns_top)
    N_lef = count_nodes(ns_left)

    # Per-node magnitudes so that Σ(node loads) = edge_length
    T_bottom = (width  / N_bot)
    T_right  = (height / N_rig)
    T_top    = (width  / N_top)
    T_left   = (height / N_lef)

    loads = LoadCollector("EDGE_TANGENTIAL")
    loads.add(CLoad(ns_bottom, (+T_bottom, 0.0     , 0.0)))  # +x
    loads.add(CLoad(ns_right,  (0.0     , -T_right, 0.0)))  # -y
    loads.add(CLoad(ns_top,    (-T_top ,  0.0     , 0.0)))  # -x
    loads.add(CLoad(ns_left,   (0.0     , +T_left , 0.0)))  # +y
    solid.add_loadcollector(loads)

    meta = {
        "width": int(width),
        "height": int(height),
        "plate_centered": True,
        "mesh_target_edge_size": 2.0,
    }
    return solid, supports, loads, hole_meta, meta


# -----------------------------------------------------------------------------
# Solution helpers
# -----------------------------------------------------------------------------
def first_buckling_value(solution) -> float:
    return float(solution.steps[0].fields["BUCKLING_FACTORS"][0].item())


def collect_hole_nodes_xy_at_midplane(model: Model) -> list[tuple[float, float]]:
    # Robust: falls kein HOLE-Nodeset existiert, leere Liste zurückgeben
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


# -----------------------------------------------------------------------------
# Single run (returns output path)
# -----------------------------------------------------------------------------
def run_single(out_dir: Path, engine_path: Path, width: int, height: int, no_holes: bool) -> Path:
    random.seed()

    model, supports, loads, hole_meta, plate_meta = build_model(width=width, height=height, no_holes=no_holes)
    model.add_step(
        LinearBucklingStep(
            num_eigenvalues=10,
            load_collectors=[loads],
            support_collectors=[supports],
        )
    )

    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
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
                "min_xy": [-plate_meta["width"] / 2.0, -plate_meta["height"] / 2.0],
                "max_xy": [ plate_meta["width"] / 2.0,  plate_meta["height"] / 2.0],
                "width": int(plate_meta["width"]),
                "height": int(plate_meta["height"]),
                "extruded_layers": 5,
                "layer_spacing": 0.2,
            },
            "material": {"E_MPa": 210e3, "nu": 0.30, "rho_t_per_mm3": 7850.0e-12},
            "loads": "outer edges tangential (+x bottom, -y right, -x top, +y left), per-node magnitude (edge_length / N_edge)",
            "supports": "BL: u1=u2=0; BR: u2=0; all outer edges: w=0",
        },
        "hole": hole_meta,
        "first_buckling_factor": float(lam1),
        "hole_nodes_xy": [[x, y] for (x, y) in hole_xy_mid],
    }
    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

    return out_path


def _run_once(seed: int, out_dir: Path, engine_path: Path, width: int, height: int, no_holes: bool):
    random.seed(seed)
    try:
        p = run_single(out_dir, engine_path, width=width, height=height, no_holes=no_holes)
        return True, seed, width, str(p)
    except Exception as e:
        return False, seed, width, str(e)


# -----------------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------------
def parse_width_seq(spec: str) -> List[int]:
    """
    Parse "start:end[:step]" into an inclusive list.
    Defaults: step=20 wenn nicht angegeben.
    Beispiel: "40:200:20" -> [40,60,80,...,200]
              "40:200"    -> [40,60,80,...,200] (step=20)
    """
    parts = spec.split(":")
    if len(parts) not in (2, 3):
        raise ValueError("width_seq must be 'start:end[:step]'")
    start = int(parts[0])
    end = int(parts[1])
    step = int(parts[2]) if len(parts) == 3 else 20
    if step <= 0:
        raise ValueError("step must be > 0")
    if end < start:
        raise ValueError("end must be >= start")
    vals = list(range(start, end + 1, step))
    if vals[-1] != end:
        vals.append(end)
    return vals


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate buckling samples with optional B-spline hole. Width can vary; height stays fixed.")
    parser.add_argument("--n", type=int, required=True, help="Number of samples to generate. If --width_seq is given: per width.")
    parser.add_argument("--out", type=Path, required=True, help="Output directory for YAML files.")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers.")
    parser.add_argument("--seed", type=int, default=None, help="Base RNG seed (optional).")
    parser.add_argument("--engine", type=Path, default=None, help="Path to FEMaster engine (overrides default).")

    # Control hole & width(s)
    parser.add_argument("--no_holes", action="store_true", help="Generate plates without a hole.")
    parser.add_argument("--width", type=int, default=40, help="Plate width for single-width runs (ignored if --width_seq is given).")
    parser.add_argument("--width_seq", type=str, default=None, help="Generate for multiple widths: 'start:end[:step]'. n is per-width.")
    parser.add_argument("--height", type=int, default=40, help="Fixed plate height (does not vary).")

    args = parser.parse_args()

    # Engine path
    default_engine = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    engine_path = args.engine if args.engine is not None else default_engine

    # Multiprocessing start method
    try:
        set_start_method("spawn")
    except RuntimeError:
        pass

    # Widths
    if args.width_seq:
        widths = parse_width_seq(args.width_seq)
    else:
        widths = [int(args.width)]

    # Seeds
    base_rng = random.Random(args.seed) if args.seed is not None else random.Random()
    total_runs = args.n * len(widths)
    seeds: List[int] = [base_rng.getrandbits(64) for _ in range(total_runs)]

    ok = fail = 0
    args.out.mkdir(parents=True, exist_ok=True)

    # Jobs aufsetzen
    jobs: List[Tuple[int, int]] = []
    for w in widths:
        for _ in range(args.n):
            jobs.append((w, 0))  # Platzhalter

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for i, (w, _) in enumerate(jobs):
            futs.append(ex.submit(_run_once, seeds[i], args.out, engine_path, w, args.height, args.no_holes))

        with tqdm.tqdm(total=len(futs), desc="Buckling runs", ncols=100, file=sys.stdout) as bar:
            for fut in as_completed(futs):
                success, seed, width, msg = fut.result()
                if success:
                    ok += 1
                else:
                    fail += 1
                    tqdm.tqdm.write(f"[seed={seed} width={width}] ERROR: {msg}", file=sys.stdout)
                bar.update(1)
                bar.set_postfix(ok=ok, fail=fail)

    print(f"\nDone. ok={ok}, failed={fail}. Output dir: {args.out.resolve()}")


if __name__ == "__main__":
    main()
