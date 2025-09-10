# sensor_opt/hbuilder.py
from __future__ import annotations
from typing import List, Dict, Tuple
import numpy as np
from .console import subhdr, kv, warn, ok
from .gauges import GAUGE_ORDER

class ChannelMap:
    """Speichert Sensor-Kanäle: (node_id, gauge_idx) -> row in H."""
    def __init__(self, pairs: List[Tuple[int, int]], labels: List[str], sigma: np.ndarray):
        self.pairs  = pairs        # List of (node_id, gauge_idx)
        self.labels = labels       # one per channel
        self.sigma  = sigma        # channel noise std (n_channels,)

def _find_strain_field_name(sol) -> str:
    """Finde ein Feld, das 'STRAIN' heißt (notfalls Case-sensitiv)."""
    lc = next(iter(sol.loadcases.keys()))
    names = list(sol.loadcases[lc].keys())
    for n in names:
        if n.upper() == "STRAIN":
            return n
    # If not present, try to guess:
    for n in names:
        if "STRAIN" in n.upper():
            warn(f"Using non-standard strain field '{n}'.")
            return n
    raise KeyError("No STRAIN field found in Solution.")

def build_channels(geometry,
                   candidate_nodes: List[int],
                   allowed_gauge_idxs: List[int],
                   sigma_by_gauge: Dict[int, float] | float = 1e-6) -> ChannelMap:
    """Erzeuge Sensor-Kanäle aus (node_id × allowed gauges)."""
    pairs: List[Tuple[int,int]] = []
    labels: List[str] = []
    if isinstance(sigma_by_gauge, (float, int)):
        sig = float(sigma_by_gauge)
        sigmas = []
    else:
        sigmas = []

    for nid in candidate_nodes:
        x, y, z = geometry.nodes[nid]
        for g in allowed_gauge_idxs:
            pairs.append((nid, g))
            labels.append(f"N{nid:06d} @({x:.3f},{y:.3f},{z:.3f}) :: {GAUGE_ORDER[g].upper()}")
            if isinstance(sigma_by_gauge, dict):
                sigmas.append(float(sigma_by_gauge.get(g, 1e-6)))

    if isinstance(sigma_by_gauge, (float, int)):
        sigma = np.full(len(pairs), sig, dtype=float)
    else:
        sigma = np.array(sigmas, dtype=float)
    return ChannelMap(pairs, labels, sigma)

def build_H_from_solution(sol,
                          geometry,
                          created_steps: List[str],
                          created_collectors: List[str],
                          channel_map: ChannelMap,
                          unitize_loads: bool = True) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """
    Baue H: Zeilen = Sensor-Kanäle, Spalten = DOF-Parameter (Reihenfolge wie created_steps).
    Jede Spalte j ist die STRAIN-Antwort des Loadcase j **pro Einheit Lastamplitude**.
    """
    subhdr("Building H from solution")
    strain_field = _find_strain_field_name(sol)

    n_ch = len(channel_map.pairs)
    n_p  = len(created_steps)
    H = np.zeros((n_ch, n_p), dtype=float)

    # Mapping Step-Reihenfolge -> loadcase key
    # Annahme: Reihenfolge der erstellten Steps == Reihenfolge in Solution (LC 1..n)
    lc_keys = sorted(sol.loadcases.keys(), key=lambda k: int(k))
    if len(lc_keys) != n_p:
        warn(f"#loadcases ({len(lc_keys)}) != #created steps ({n_p}). Using min().")
    n_use = min(len(lc_keys), n_p)

    for j in range(n_use):
        lc = lc_keys[j]
        S = sol.get(lc, strain_field)  # shape (n_nodes, 6)
        # falls bei dir elementar: hier ggf. Element->Knoten-Map einziehen

        # Amplitude der einzigen DOF-Komponente dieses Steps
        # -> hole sie aus geometry.loads[collector]['data']
        collector = created_collectors[j]
        data = geometry.loads[collector]["data"]
        amp = next((v for v in data if v is not None), 1.0)
        if unitize_loads and amp not in (0, None):
            scale = 1.0 / float(amp)
        else:
            scale = 1.0

        # befülle Spalte j
        for r, (nid, g) in enumerate(channel_map.pairs):
            # Annahme: row index == node_id
            if nid >= S.shape[0]:
                raise IndexError(f"Node id {nid} out of range for STRAIN field with shape {S.shape}.")
            H[r, j] = S[nid, g] * scale

    ok(f"H built: shape = {H.shape}")
    return H, channel_map.labels, channel_map.sigma
