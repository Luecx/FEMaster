# file: run_buckle_with_circular_bspline_hole.py
from __future__ import annotations

"""
Buckling of a 40×40 plate (about the origin) with a circular-ish B-spline hole.

This reproduces the Abaqus-style randomness and sizing from your CAE snippet:
- Plate outline: rectangle from (-20, -20) to (20, 20).
- Hole control points:
    n_ctrl  ~ randint(3, 5)
    radius  ~ uniform(4, 12)
    per-point jitter in x and y: uniform(-2, 2)
- Target in-plane edge size ≈ 2.0 → 20 subdivisions per side.

Modeling:
- 2D quad-ish mesh with one inner loop (the hole) → extruded to thin solid:
  thickness = 1.0 (5 layers × 0.2).
- Supports:
    BL corner: u1=0, u2=0
    BR corner:       u2=0
    All outer edges: w=0
- Loads (tangential, magnitude 1.0 on outer edges):
    bottom:+x, right:-y, top:-x, left:+y

Outputs:
- YAML file: buckling_bspline_run_YYYYMMDD-HHMMSS.yaml
  Contains:
    * meta, hole (radius + control points),
    * first_buckling_factor,
    * hole_nodes_xy: unique (x, y) of HOLE nodes at z = 0 (no duplicates).
"""

import math
import random
import time
import tqdm
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import set_start_method

if __package__ in (None, ""):
    import os
    import sys
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

import yaml


# -----------------------------------------------------------------------------
# Geometry helpers (Abaqus-style randomness)
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
    """
    Create a closed B-spline approximating a circle, mirroring the CAE logic.

    Returns:
        bspline, control_points[(x, y)], radius
    """
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
def build_model_with_bspline_hole() -> tuple[Model, SupportCollector, LoadCollector, dict]:
    """
    Build thin solid plate with a HOLE B-spline, quad-ish mesh, and buckling BCs.
    """
    # --- Outer rectangle (20 subdivisions per side → ~2.0 size over 40 length)
    seg_bottom = StraightSegment((-20.0, -20.0), (+20.0, -20.0), n_subdivisions=20, name="BOTTOM")
    seg_right  = StraightSegment((+20.0, -20.0), (+20.0, +20.0), n_subdivisions=20, name="RIGHT")
    seg_top    = StraightSegment((+20.0, +20.0), (-20.0, +20.0), n_subdivisions=20, name="TOP")
    seg_left   = StraightSegment((-20.0, +20.0), (-20.0, -20.0), n_subdivisions=20, name="LEFT")
    outer = SegmentGroup([seg_bottom, seg_right, seg_top, seg_left], name="OUTER")

    # --- Hole: Abaqus-style random B-spline (center at origin by construction)
    bspline, ctrl_pts, radius = make_abaqus_like_bspline(name="HOLE")
    hole = SegmentGroup([bspline], name="HOLE")

    # --- Mesh (2D → quad-ish), then extrude to thin solid of thickness 1.0
    plate2d = Model.mesh_2d([outer, hole], mesh_type=2, name="Plate40x40_withBSplineHole_2D")
    solid = plate2d.extruded(n=5, spacing=0.2)  # 5 layers × 0.2 → total thickness = 1.0
    solid.name = "SolidPlate_40x40_t1p0_n5_bspline"

    # --- Elements & section
    all_elems = [e for e in solid.elements._items if e is not None]
    elem_set = ElementSet("EALL", all_elems)
    solid.add_elementset(elem_set)

    steel = solid.add_material(
        Material("STEEL")
        .set_elasticity(ElasticityIsotropic(210e3, 0.30))
        .set_density(Density(7850.0e-12))
    )
    solid.add_section(SolidSection(material=steel, elset=elem_set))

    # --- Node sets (edges auto-named; create persistent corner sets)
    ns_bottom = solid.node_sets["BOTTOM"]
    ns_right  = solid.node_sets["RIGHT"]
    ns_top    = solid.node_sets["TOP"]
    ns_left   = solid.node_sets["LEFT"]

    ns_bl = ns_bottom & ns_left   # bottom-left intersection
    ns_br = ns_bottom & ns_right  # bottom-right intersection
    ns_bl.name = "BOTTOM_LEFT"
    ns_br.name = "BOTTOM_RIGHT"
    solid.add_nodeset(ns_bl)
    solid.add_nodeset(ns_br)

    # --- Supports
    supports = SupportCollector("SUPPORTS")
    supports.add(Support(ns_bl, (0.0, 0.0, None, None, None, None)))  # BL: u1=0, u2=0
    supports.add(Support(ns_br, (None, 0.0, None, None, None, None))) # BR:       u2=0
    for edge in (ns_bottom, ns_right, ns_top, ns_left):
        supports.add(Support(edge, (None, None, 0.0, None, None, None)))  # w=0 on edges
    solid.add_supportcollector(supports)

    # --- Tangential edge loads ±1.0 on opposite edges
    loads = LoadCollector("EDGE_TANGENTIAL")
    T = 1.0 / 6.0  # magnitude
    loads.add(CLoad(ns_bottom, (+  T,  0.0, 0.0)))  # +x
    loads.add(CLoad(ns_right,  ( 0.0, -  T, 0.0)))  # -y
    loads.add(CLoad(ns_top,    (-  T,  0.0, 0.0)))  # -x
    loads.add(CLoad(ns_left,   ( 0.0, +  T, 0.0)))  # +y
    solid.add_loadcollector(loads)

    # --- Metadata for output
    hole_meta = {
        "radius": float(radius),
        "control_points_xy": [[float(x), float(y)] for (x, y) in ctrl_pts],
        "note": "Abaqus-like: uniform angles, per-point ±2 x/y jitter",
    }
    return solid, supports, loads, hole_meta


# -----------------------------------------------------------------------------
# Solution helpers (no fallbacks, no generalization)
# -----------------------------------------------------------------------------
def first_buckling_value(solution) -> float:
    """Return λ₁ from the official location: steps[0].fields['EIGENVALUES'][0]."""
    return float(solution.steps[0].fields["BUCKLING_FACTORS"][0].item())


def collect_hole_nodes_xy_at_midplane(model: Model) -> list[tuple[float, float]]:
    """
    Collect unique (x, y) of HOLE nodes restricted to the z = 0 layer.
    - No general cases, no None-guards.
    - Uniqueness enforced on (x, y) with rounding.
    """
    ns_hole = model.node_sets["HOLE"]

    out_xy: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()

    for n in ns_hole:  # NodeSet is directly iterable over Node objects
        # node has attributes: n.x, n.y, n.z (and n.node_id)
        if abs(float(n.z) - 0.0) < 1e-12:
            xy = (round(float(n.x), 12), round(float(n.y), 12))
            if xy not in seen:
                seen.add(xy)
                out_xy.append((float(n.x), float(n.y)))

    return out_xy


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    # Keep variability like the Abaqus snippet
    random.seed()

    model, supports, loads, hole_meta = build_model_with_bspline_hole()

    # Buckling step
    model.add_step(
        LinearBucklingStep(
            num_eigenvalues=10,
            load_collectors=[loads],
            support_collectors=[supports],
        )
    )

    engine_path = Path(__file__).resolve().parent.parent / "bin" / "FEMaster"
    runner = Runner().set_model(model)
    runner.set_engine(Runner.Engine.FEMASTER, path=engine_path)
    solution = runner.run()

    lam1 = first_buckling_value(solution)
    hole_xy_mid = collect_hole_nodes_xy_at_midplane(model)

    # Dump YAML (flat)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    out_path = Path(f"buckling_bspline_run_{timestamp}.yaml")
    payload = {
        "meta": {
            "timestamp": timestamp,
            "plate": {
                "min_xy": [-20.0, -20.0],
                "max_xy": [ 20.0,  20.0],
                "extruded_layers": 5,
                "layer_spacing": 0.2,
            },
            "material": {"E": 210e3, "nu": 0.30, "rho": 7850.0},
            "loads": "outer edges tangential (+x bottom, -y right, -x top, +y left), magnitude 1.0",
            "supports": "BL: u1=u2=0; BR: u2=0; all outer edges: w=0",
        },
        "hole": hole_meta,
        "first_buckling_factor": float(lam1),
        "hole_nodes_xy": [[x, y] for (x, y) in hole_xy_mid],  # z=0 only, xy unique
    }
    with open(out_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)


def _run_once(seed: int):
    random.seed(seed)
    try:
        main()
        return True, seed, ""
    except Exception as e:
        return False, seed, str(e)

def run_all():
    N = 10000
    WORKERS = 8
    base_rng = random.Random()
    seeds = [base_rng.getrandbits(64) for _ in range(N)]

    ok = fail = 0
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(_run_once, s) for s in seeds]
        with tqdm.tqdm(total=N, desc="Buckling runs", ncols=100, file=sys.stdout) as bar:
            for fut in as_completed(futs):
                success, seed, msg = fut.result()   # now matches
                if success:
                    ok += 1
                else:
                    fail += 1
                    tqdm.tqdm.write(f"seed={seed}: {msg}", file=sys.stdout)
                bar.update(1)
                bar.set_postfix(ok=ok, fail=fail)

    print(f"\nDone. ok={ok}, failed={fail}")

if __name__ == "__main__":
    run_all()