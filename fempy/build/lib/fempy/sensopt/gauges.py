# sensor_opt/gauges.py
from __future__ import annotations
from typing import Dict, List

# Standardreihenfolge im STRAIN-Feld: exx, eyy, ezz, eyz, ezx, exy
GAUGE_ORDER = ["exx", "eyy", "ezz", "eyz", "ezx", "exy"]

# Synonyme → Index
GAUGE_TO_INDEX: Dict[str, int] = {
    "exx": 0, "strain_xx": 0, "eps_xx": 0, "εxx": 0, "xx": 0, "x": 0,
    "eyy": 1, "strain_yy": 1, "eps_yy": 1, "εyy": 1, "yy": 1, "y": 1,
    "ezz": 2, "strain_zz": 2, "eps_zz": 2, "εzz": 2, "zz": 2, "z": 2,
    "eyz": 3, "strain_yz": 3, "eps_yz": 3, "εyz": 3, "yz": 3,
    "ezx": 4, "strain_zx": 4, "eps_zx": 4, "εzx": 4, "zx": 4,
    "exy": 5, "strain_xy": 5, "eps_xy": 5, "εxy": 5, "xy": 5,
}

def parse_allowed_gauges(names: List[str]) -> List[int]:
    """Wandle eine Liste von Namen/Synonymen in eindeutige Gauge-Indizes (0..5) um."""
    idxs = []
    for n in names:
        key = str(n).strip().lower()
        if key not in GAUGE_TO_INDEX:
            raise ValueError(f"Unknown gauge name '{n}'. Allowed: {GAUGE_ORDER}")
        idx = GAUGE_TO_INDEX[key]
        if idx not in idxs:
            idxs.append(idx)
    return idxs
