# sensor_opt/select.py
from __future__ import annotations
from typing import Optional, List, Dict
import numpy as np
from .console import subhdr, kv, ok
from .fim_selector import SelectorConfig, SensorSelector  # <— hier angepasst!

def run_selection(H: np.ndarray,
                  sigma: np.ndarray | float,
                  param_std_targets: List[float],
                  labels: List[str],
                  removal_criterion: str = "D",
                  cond_max: float = 1e12,
                  prior_prec: float = 1e-12):
    """
    Führt die FIM-basierte Backward-Elimination aus und gibt das Ergebnis zurück.
    """
    subhdr("FIM-based sensor selection")

    if isinstance(sigma, (float, int)):
        sigma_vec = np.full(H.shape[0], float(sigma))
    else:
        sigma_vec = np.asarray(sigma, float)

    cfg = SelectorConfig(
        H=H,
        sigma=sigma_vec,
        param_std_max=np.asarray(param_std_targets, float),
        removal_criterion=removal_criterion,
        cond_max=cond_max,
        prior_prec=prior_prec,
        sensor_labels=labels,
    )
    sel = SensorSelector(cfg)
    result = sel.run()
    ok("Selection finished.")
    return result
