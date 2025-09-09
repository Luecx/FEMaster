# sensor_opt/runner.py
from __future__ import annotations
import os, subprocess, shutil
from typing import Tuple
from .console import hdr, kv, ok, err
from ..solution import Solution  # deine Solution-Klasse

class SensorRunner:
    def __init__(self, geometry, solver_path: str, work_dir: str,
                 device: str = "gpu", method: str = "indirect", ncpus: int = 1):
        self.geometry  = geometry
        self.solver    = solver_path
        self.work_dir  = work_dir
        self.device    = device
        self.method    = method
        self.ncpus     = ncpus

    def _write_input(self, fname: str) -> str:
        os.makedirs(self.work_dir, exist_ok=True)
        path = os.path.join(self.work_dir, fname)
        self.geometry.write_input_deck(path)
        return path

    def run(self, model_name: str = "sensoropt_model") -> Tuple[str, Solution]:
        hdr("Running FE solver for SensorOpt")
        inppath = self._write_input(model_name + ".inp")
        kv("Input deck", inppath)
        log_path = os.path.join(self.work_dir, model_name + ".log")
        ok("Launching solver ...")

        with open(log_path, 'w') as log:
            res = subprocess.run([self.solver, inppath, "--ncpus", str(self.ncpus)],
                                 stdout=log, stderr=log)

        if res.returncode != 0:
            with open(log_path, 'r') as f:
                logtxt = f.read()[-2000:]
            err("Solver returned non-zero exit code.")
            raise RuntimeError(f"Solver failed. Tail log:\n{logtxt}")

        res_path = os.path.join(self.work_dir, model_name + ".res")
        if not os.path.exists(res_path):
            raise FileNotFoundError(f"Result file not found: {res_path}")

        ok("Reading solution file ...")
        sol = Solution.open(res_path, loadingbar=False)
        kv("Loadcases", len(sol.loadcases))
        return res_path, sol
