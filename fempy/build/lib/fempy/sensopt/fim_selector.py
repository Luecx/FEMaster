# sensor_opt/fim_selector.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Iterable, Optional, Tuple, List
import numpy as np

# -------------------- Helfer --------------------

def _as_sym_posdef(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, float)
    M = 0.5 * (M + M.T)
    try:
        np.linalg.cholesky(M)
    except np.linalg.LinAlgError:
        eps = 1e-12 * max(1.0, np.linalg.norm(M))
        M = M + eps * np.eye(M.shape[0])
    return M

def _pretty_bar(value: float, limit: float, width: int = 28, invert: bool = False) -> str:
    if limit <= 0: limit = 1.0
    t = value / limit
    t = max(0.0, min(1.0, t))
    if invert: t = 1.0 - t
    n = int(round(t * width))
    return "[" + "#" * n + "-" * (width - n) + "]"

# -------------------- Config --------------------

@dataclass
class SelectorConfig:
    # Designmatrix / Sensitivitäten (n_sensors x p)
    H: Optional[np.ndarray] = None

    # Rauschen (genau eine Angabe)
    sigma: Optional[np.ndarray] = None     # (n,)  diag-Σ
    Sigma: Optional[np.ndarray] = None     # (n,n) voll

    # Ziele
    param_std_max: Optional[np.ndarray] = None  # (p,)
    qoi_gradients: Optional[np.ndarray] = None  # (n_qoi, p)
    qoi_var_max: Optional[float] = None

    # Algorithmus
    removal_criterion: str = "D"   # "D" oder "A"
    prior_prec: float = 1e-12
    cond_max: float = 1e12

    # Ausgabe
    sensor_labels: Optional[Iterable] = None

# -------------------- Selektor --------------------

class SensorSelector:
    """
    Backward-Elimination anhand Fisher-Informationsmatrix:
      - D- oder A-Optimalität
      - diagonales oder volles Messrauschen
      - feste Grenzen auf std(theta), Konditionszahl(F), optional Var[QoI]
    """

    def __init__(self, config: SelectorConfig):
        self.cfg = config
        self._validate_and_prepare()

    def _validate_and_prepare(self):
        cfg = self.cfg
        if cfg.H is None:
            raise ValueError("SelectorConfig.H muss gesetzt sein.")
        H = np.asarray(cfg.H, float)
        if H.ndim != 2:
            raise ValueError("H muss 2D sein (n_sensors x p).")
        n, p = H.shape

        # Rauschen
        if (cfg.sigma is None) == (cfg.Sigma is None):
            raise ValueError("Entweder sigma (diag) ODER Sigma (voll) angeben.")
        if cfg.sigma is not None:
            sigma = np.asarray(cfg.sigma, float).reshape(-1)
            if sigma.shape[0] != n:
                raise ValueError("sigma hat falsche Länge.")
            self._diag_sigma = True
            self._sigma = sigma
            Sigma_inv = np.diag(1.0 / (sigma**2))
        else:
            Sigma = np.asarray(cfg.Sigma, float)
            if Sigma.shape != (n, n):
                raise ValueError("Sigma muss (n x n) sein.")
            Sigma = _as_sym_posdef(Sigma)
            Sigma_inv = np.linalg.inv(Sigma)
            # falls faktisch diagonal, aktiviere Schnellpfad
            self._diag_sigma = np.allclose(Sigma, np.diag(np.diag(Sigma)))
            if self._diag_sigma:
                self._sigma = np.sqrt(np.diag(Sigma))

        # Ziele
        if cfg.param_std_max is None:
            raise ValueError("param_std_max muss gesetzt sein.")
        param_std_max = np.asarray(cfg.param_std_max, float).reshape(-1)
        if param_std_max.shape[0] != p:
            raise ValueError("param_std_max falsche Länge.")

        # QoI
        G = None
        if cfg.qoi_gradients is not None and cfg.qoi_var_max is not None:
            G = np.asarray(cfg.qoi_gradients, float)
            if G.ndim != 2 or G.shape[1] != p:
                raise ValueError("qoi_gradients muss (n_qoi x p) sein.")

        # Labels
        if cfg.sensor_labels is None:
            labels = [f"idx={i}" for i in range(n)]
        else:
            labels = [str(x) for x in cfg.sensor_labels]
            if len(labels) != n:
                raise ValueError("sensor_labels falsche Länge.")

        # Kriterium
        crit = cfg.removal_criterion.upper()
        if crit not in ("D", "A"):
            raise ValueError("removal_criterion muss 'D' oder 'A' sein.")

        self.H = H
        self.Sigma_inv = Sigma_inv
        self.param_std_max = param_std_max
        self.G = G
        self.crit = crit
        self.p = p
        self.n = n
        self.labels = labels
        self.prior = float(cfg.prior_prec)
        self.cond_max = float(cfg.cond_max)

    # ---------- Mathe ----------

    def _fim_from_active(self, active: List[int]) -> np.ndarray:
        idx = np.array(active, dtype=int)
        H_a = self.H[idx, :]
        Sii = self.Sigma_inv[np.ix_(idx, idx)]
        F = H_a.T @ Sii @ H_a
        if self.prior:
            F = F + self.prior * np.eye(self.p)
        return _as_sym_posdef(F)

    def _u_rank1(self, j: int) -> np.ndarray:
        # Beitrag dieses Sensors bei diagonalem Σ: u = H_j / sigma_j
        return (self.H[j, :] / float(self._sigma[j])).reshape(-1, 1)

    def _choose_removals(self, F, Finv, active, k):
        """Wählt die k 'besten' Sensoren zum Entfernen."""
        if self._diag_sigma:
            costs = []
            for j in active:
                u = self._u_rank1(j)
                s = (u.T @ Finv @ u).item()
                if 1.0 - s <= 1e-12:
                    val = np.inf
                else:
                    if self.crit == "D":
                        val = -np.log(1.0 - s)
                    else:
                        Finv_u = Finv @ u
                        val = ((Finv_u.T @ Finv_u) / (1.0 - s)).item()
                costs.append((val, j))
            # sortiere nach kleinstem Verlust
            costs.sort(key=lambda x: x[0])
            return [j for _, j in costs[:k]]
        else:
            base_logdet = np.linalg.slogdet(F)[1] if self.crit == "D" else None
            base_trace  = np.trace(np.linalg.inv(F)) if self.crit == "A" else None
            costs = []
            aset = set(active)
            for j in active:
                idx_next = sorted(list(aset - {j}))
                F_new = self._fim_from_active(idx_next)
                if self.crit == "D":
                    drop = (base_logdet - np.linalg.slogdet(F_new)[1]).item()
                    val = max(drop, 0.0)
                else:
                    inc = (np.trace(np.linalg.inv(F_new)) - base_trace).item()
                    val = max(inc, 0.0)
                costs.append((val, j))
            costs.sort(key=lambda x: x[0])
            return [j for _, j in costs[:k]]


    def _downdate_inverse(self, F: np.ndarray, Finv: np.ndarray, u: np.ndarray) -> Tuple[np.ndarray, np.ndarray, bool]:
        s = (u.T @ Finv @ u).item()
        if 1.0 - s <= 1e-12:
            return F, Finv, False
        Finv_u = Finv @ u
        Finv_new = Finv + (Finv_u @ Finv_u.T) / (1.0 - s)
        F_new = F - (u @ u.T)
        return _as_sym_posdef(F_new), _as_sym_posdef(Finv_new), True

    # ---------- Ausgabe ----------

    def _print_header(self, active: List[int]):
        print("\n" + "=" * 82)
        print("BLACKBOX SENSOR-ELIMINATION (Backward)")
        print("=" * 82)
        print(f"Parameter p={self.p} | Start-Sensoren n={self.n} | Kriterium: {self.crit}-optimal")
        print(f"Std-Grenzen: {', '.join(f'{s:.4g}' for s in self.param_std_max)}")
        if self.G is not None and self.cfg.qoi_var_max is not None:
            print(f"QoI-Varianz-Grenze: <= {self.cfg.qoi_var_max:.4g}")
        print(f"Konditionszahl-Grenze: <= {self.cond_max:.1e}")
        print("-" * 82)
        print("Sensor-Kandidaten:")
        print("\n".join(self.labels[i] for i in active))
        print("-" * 82)

    def _print_step(self, step: int, j: int, cost: float, metric_name: str,
                    n_active_next: int, stds: np.ndarray, condF: float, max_qoi_var: Optional[float]):
        print(f"\nStep {step:02d} | entferne {self.labels[j]} | {metric_name} = {cost:.6g}")
        print(f"  aktive Sensoren: {n_active_next}")
        print(f"  std(theta):      " + ", ".join(f"{s:.4g}" for s in stds))
        print(f"  Konditionszahl:  {condF:.3e}  {_pretty_bar(min(condF,self.cond_max), self.cond_max)} (<= {self.cond_max:.1e})")
        if max_qoi_var is not None:
            print(f"  max Var[QoI]:    {max_qoi_var:.4g}")

    # ---------- Hauptlauf ----------

    def run(self, start_percent=0.1):
        H = self.H
        n, p = H.shape
        active = list(range(n))
        F = self._fim_from_active(active)
        Finv = np.linalg.inv(F)

        self._print_header(active)
        step = 0

        while len(active) > p:
            # Initiale Blockgröße in Sensoren
            block_size = max(1, int(len(active) * start_percent))

            removed_this_round = False
            while block_size >= 1:
                to_remove = self._choose_removals(F, Finv, active, block_size)
                next_active = [i for i in active if i not in to_remove]
                F_new = self._fim_from_active(next_active)
                Finv_new = np.linalg.inv(F_new)

                stds  = np.sqrt(np.diag(Finv_new))
                condF = float(np.linalg.cond(F_new))
                ok_params = np.all(stds <= self.param_std_max + 1e-12)
                ok_cond   = (condF <= self.cond_max)

                ok_qoi, max_qoi_var = True, None
                if self.G is not None and self.cfg.qoi_var_max is not None:
                    S = self.G @ Finv_new @ self.G.T
                    max_qoi_var = float(np.max(np.diag(S)))
                    ok_qoi = (max_qoi_var <= self.cfg.qoi_var_max + 1e-12)

                if ok_params and ok_cond and ok_qoi:
                    step += 1
                    names = [self.labels[j] for j in to_remove]
                    if len(names) > 4:
                        display_names = names[:3] + ["..."] + [names[-1]]
                    else:
                        display_names = names

                    print(f"\nStep {step:02d} | entferne {block_size} Sensor(en): "
                          + ", ".join(display_names))

                    print(f"  aktive Sensoren: {len(next_active)}")
                    print(f"  std(theta):      " + ", ".join(f"{s:.4g}" for s in stds))
                    print(f"  Konditionszahl:  {condF:.3e}  "
                          f"{_pretty_bar(min(condF,self.cond_max), self.cond_max)} (<= {self.cond_max:.1e})")
                    if max_qoi_var is not None:
                        print(f"  max Var[QoI]:    {max_qoi_var:.4g}")

                    # Übernehmen
                    active = next_active
                    F, Finv = F_new, Finv_new
                    removed_this_round = True
                    break
                else:
                    block_size //= 2  # Block verkleinern und erneut versuchen

            if not removed_this_round:
                print("\n-- STOP: Keine weitere Entfernung möglich ohne Ziele zu verletzen. --")
                break

        # Abschlussausgabe
        stds  = np.sqrt(np.diag(Finv))
        condF = float(np.linalg.cond(F))
        chosen_labels = [self.labels[i] for i in active]
        print("\n" + "-" * 82)
        print("GEWÄHLTE SENSOREN:")
        print("\n".join(f"- {lab}" for lab in chosen_labels))
        print("-" * 82)
        print("Finale std(theta): " + ", ".join(f"{s:.6g}" for s in stds))
        print(f"Finale Konditionszahl(F): {condF:.3e}")
        print("=" * 82 + "\n")

        return {
            "active_indices": active,
            "active_labels": chosen_labels,
            "F": F,
            "Cov": Finv,
            "std_theta": stds,
            "condF": condF,
        }

