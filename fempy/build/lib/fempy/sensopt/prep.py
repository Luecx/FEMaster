# sensor_opt/prep.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Iterable, Tuple, List, Dict, Optional
import copy, os
import numpy as np

from .console import hdr, subhdr, kv, bullet, ok, warn
from .gauges import parse_allowed_gauges, GAUGE_ORDER
from .runner import SensorRunner
from .hbuilder import build_channels, build_H_from_solution, ChannelMap
from .select import run_selection

@dataclass
class SensorOptConfig:
    solver_path: str
    work_dir: str = "./sensoropt_out"
    device: str = "cpu"
    method: str = "direct"
    ncpus: int = 1

    candidate_node_set: str = "NALL"
    surface_only: bool = True
    allowed_gauges: List[str] = field(default_factory=lambda: ["exx", "eyy", "ezz", "eyz", "ezx", "exy"])
    sigma_by_gauge: float | Dict[int, float] = 1e-6  # pro Kanalrauschen

    # save selected sensor nodes here
    selected_node_set: str = "SENSOR_NODES"

    # Selektion (optional)
    removal_criterion: str = "D"
    cond_max: float = 1e12
    prior_prec: float = 1e-12

class SensorOpt:
    """
    End-to-end Hilfsklasse:
      1) Steps aufbereiten (DOF-Splitting je beobachtetem Load-Collector)
      2) Kandidaten (Oberfläche/allgemein) + Kanäle (Gauges) festlegen
      3) Solver ausführen, Solution lesen
      4) H-Matrix bauen (Unit-Responses)
      5) (optional) Selektion per FIM
    """
    DOF_NAMES = ("FX", "FY", "FZ", "MX", "MY", "MZ")

    def __init__(self,
                 geometry,
                 watch_load_collectors: Iterable[Tuple[str, float]],
                 support_collectors: Optional[Iterable[str]],
                 config: SensorOptConfig,
                 in_place: bool = True):
        self.orig_geometry = geometry
        self.geometry = geometry if in_place else copy.deepcopy(geometry)
        self.watch = list(watch_load_collectors)  # [(name, allowed_sigma_param), ...]
        self.support_collectors = support_collectors
        self.cfg = config

        # Preparation outputs
        self.candidate_nodes: List[int] = []
        self.allowed_gauge_idxs: List[int] = parse_allowed_gauges(self.cfg.allowed_gauges)
        self.channel_map: Optional[ChannelMap] = None

        self.created_collectors: List[str] = []
        self.created_steps: List[str] = []
        self.split_map: Dict[str, List[str]] = {}
        self.param_targets: List[float] = []  # std targets in Spaltenreihenfolge (created steps)

        # Post-solve
        self.H: Optional[np.ndarray] = None
        self.labels: Optional[List[str]] = None
        self.sigma_channels: Optional[np.ndarray] = None

    # ---------- Pipeline ----------

    def prepare(self) -> "SensorOpt":
        hdr("SensorOpt: preparation")
        self._choose_support()
        self._compute_candidates()
        self._clear_steps()
        self._split_collectors()
        self._build_channels()
        self._print_summary()
        return self

    def solve(self, model_name: str = "sensoropt_model"):
        runner = SensorRunner(self.geometry,
                              solver_path=self.cfg.solver_path,
                              work_dir=self.cfg.work_dir,
                              device=self.cfg.device,
                              method=self.cfg.method,
                              ncpus=self.cfg.ncpus)
        res_path, sol = runner.run(model_name=model_name)
        return res_path, sol

    def build_H(self, sol):
        H, labels, sigma = build_H_from_solution(
            sol=sol,
            geometry=self.geometry,
            created_steps=self.created_steps,
            created_collectors=self.created_collectors,
            channel_map=self.channel_map,
            unitize_loads=True
        )
        self.H, self.labels, self.sigma_channels = H, labels, sigma
        return H, labels, sigma

    # --- inside class SensorOpt ---
    def select(self):
        if self.H is None:
            raise RuntimeError("Call build_H() before select().")
        result = run_selection(
            H=self.H,
            sigma=self.sigma_channels,
            param_std_targets=self.param_targets,
            labels=self.labels,
            removal_criterion=self.cfg.removal_criterion,
            cond_max=self.cfg.cond_max,
            prior_prec=self.cfg.prior_prec,
        )

        # NEW: persist selected node IDs into a node set
        self._save_selected_nodes(result["active_indices"])
        return result

    # NEW helper
    def _save_selected_nodes(self, active_indices: List[int]):
        if self.channel_map is None:
            warn("Channel map missing; cannot save selected nodes.")
            return
        # active_indices refer to H rows ⇒ map back to node IDs via channel_map.pairs
        node_ids = sorted({ self.channel_map.pairs[i][0] for i in active_indices })

        setname = self.cfg.selected_node_set
        # overwrite or create the set
        if setname in self.geometry.node_sets:
            self.geometry.node_sets[setname] = []
        else:
            self.geometry.add_node_set(setname)

        for nid in node_ids:
            self.geometry.add_node_to_set(setname, nid)

        kv("Saved selected sensor nodes", f"{setname} ({len(node_ids)})")


    # ---------- Internals ----------

    def _choose_support(self):
        if self.support_collectors:
            for sc in self.support_collectors:
                if sc not in self.geometry.supps:
                    raise ValueError(f"Support collector '{sc}' not found.")
                kv("Support collector", sc)
            return
        if len(self.geometry.supps) == 1:
            self.support_collector = next(iter(self.geometry.supps.keys()))
            kv("Support collector (auto)", self.support_collector)
            return
        raise ValueError("Ambiguous supports. Please specify 'support_collector'.")

    def _compute_candidates(self):
        # Knoten aus gewünschtem Set
        set_nodes = set(self.geometry.node_sets.get(self.cfg.candidate_node_set, []))
        if self.cfg.surface_only:
            # reine lonely-face-Oberfläche aus Geometry
            surf_nodes = set(self.geometry.get_surface_nodes())
            self.candidate_nodes = sorted(set_nodes & surf_nodes)
            kv("Candidates", f"{len(self.candidate_nodes)} nodes (intersection of '{self.cfg.candidate_node_set}' with surface)")
        else:
            self.candidate_nodes = sorted(set_nodes)
            kv("Candidates", f"{len(self.candidate_nodes)} nodes (set='{self.cfg.candidate_node_set}')")


    def _clear_steps(self):
        n_old = len(self.geometry.steps)
        self.geometry.steps = []
        kv("Cleared steps", n_old)

    def _unique_name(self, base: str, existing: set) -> str:
        if base not in existing:
            return base
        k = 1
        while f"{base}_{k}" in existing:
            k += 1
        return f"{base}_{k}"

    def _split_collectors(self):
        subhdr("Splitting load collectors into DOF components")
        existing = set(self.geometry.loads.keys())
        self.created_collectors.clear()
        self.created_steps.clear()
        self.split_map.clear()
        self.param_targets.clear()

        for orig_name, allowed_std in self.watch:
            if orig_name not in self.geometry.loads:
                raise ValueError(f"Load collector '{orig_name}' not found.")
            entry = self.geometry.loads[orig_name]
            set_name = entry["set"]
            comps = list(entry["data"])  # (fx, fy, fz, mx, my, mz)

            created_for_this: List[str] = []
            for idx, val in enumerate(comps):
                if val is None:
                    continue
                dof = self.DOF_NAMES[idx]
                new_lc = self._unique_name(f"SENSOPT_{orig_name}_{dof}", existing)
                existing.add(new_lc)

                only = [None]*6
                only[idx] = val
                self.geometry.add_load(new_lc, set_name,
                                       fx=only[0], fy=only[1], fz=only[2],
                                       mx=only[3], my=only[4], mz=only[5])
                self.created_collectors.append(new_lc)
                created_for_this.append(new_lc)

                step_name = f"SENSOPT-{orig_name}-{dof}"
                self.geometry.add_step(type="LINEAR STATIC",
                                       name=step_name,
                                       loads=[new_lc],
                                       supps=self.support_collectors)
                self.created_steps.append(step_name)
                self.param_targets.append(float(allowed_std))  # gleicher Grenzwert je DOF des Collectors

            self.split_map[orig_name] = created_for_this

        kv("Created collectors", len(self.created_collectors))
        kv("Created steps", len(self.created_steps))

    def _build_channels(self):
        subhdr("Building sensor channels from candidates & gauges")
        idxs = self.allowed_gauge_idxs
        gauge_list = [GAUGE_ORDER[i] for i in idxs]
        kv("Allowed gauges", ", ".join(gauge_list))
        self.channel_map = build_channels(
            geometry=self.geometry,
            candidate_nodes=self.candidate_nodes,
            allowed_gauge_idxs=idxs,
            sigma_by_gauge=self.cfg.sigma_by_gauge
        )
        kv("Channels", f"{len(self.channel_map.pairs)} (nodes × gauges)")

    def _print_summary(self):
        subhdr("Preparation summary")
        kv("Work dir", os.path.abspath(self.cfg.work_dir))
        kv("Solver", self.cfg.solver_path)
        kv("Steps total", len(self.created_steps))
        kv("Param targets", f"{len(self.param_targets)} (std limits)")
        bullet([f"{k} -> {len(v)} DOFs" for k, v in self.split_map.items()])
        ok("Preparation finished.")
