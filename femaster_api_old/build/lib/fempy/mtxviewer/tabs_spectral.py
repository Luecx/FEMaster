# =============================
# File: mtxviewer/tabs_spectral.py
# =============================
from typing import Optional, Tuple
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, svds
from PyQt5 import QtCore, QtWidgets


class SpectralTab(QtWidgets.QWidget):
    """
    Spectral analysis on the symmetric part S = (A + A^T)/2 (float64 CSC).

    - Smallest algebraic eigenvalues (LOW): robustly estimate λ_min, then
      use shift-invert around that estimate to retrieve the k smallest.
    - Largest algebraic eigenvalues (HIGH): eigsh(S, which='LA') -> sort desc.
    - Top singular values (SVD): svds(A) -> sort desc.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self.worker: Optional[QtCore.QThread] = None
        self._build()

    def _build(self):
        v = QtWidgets.QVBoxLayout(self)
        ctl = QtWidgets.QHBoxLayout()
        self.k_spin = QtWidgets.QSpinBox(); self.k_spin.setRange(1, 500); self.k_spin.setValue(10)
        self.btn_low  = QtWidgets.QPushButton("Smallest eigenvalues (ascending)");  self.btn_low.setProperty('primary', True)
        self.btn_high = QtWidgets.QPushButton("Largest eigenvalues (descending)");   self.btn_high.setProperty('primary', True)
        # self.btn_svd  = QtWidgets.QPushButton("Top singular values (descending)");   self.btn_svd.setProperty('primary', True)
        ctl.addWidget(QtWidgets.QLabel("k:")); ctl.addWidget(self.k_spin); ctl.addStretch(1)
        ctl.addWidget(self.btn_low);
        ctl.addWidget(self.btn_high);
        #ctl.addWidget(self.btn_svd)
        v.addLayout(ctl)

        self.table = QtWidgets.QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["index", "value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table)

        self.log = QtWidgets.QPlainTextEdit(); self.log.setReadOnly(True); self.log.setMaximumHeight(120)
        v.addWidget(self.log)

        self.btn_low.clicked.connect(lambda: self._compute('LOW'))
        self.btn_high.clicked.connect(lambda: self._compute('HIGH'))
        # self.btn_svd.clicked.connect(lambda: self._compute('SVD'))

    def set_matrix(self, A: sparse.spmatrix):
        self.A = A
        self.table.setRowCount(0)
        self.log.setPlainText("")

    def _compute(self, which: str):
        if self.A is None: return
        k = int(self.k_spin.value())
        A = self.A

        def task():
            if which in ('LOW', 'HIGH'):
                # Work on symmetric part
                S = ((A + A.T) * 0.5).astype(np.float64).tocsc()
                n = S.shape[0]
                kk = int(min(k, max(1, n - 2)))

                if which == 'HIGH':
                    vals = eigsh(S, k=kk, which='LA', return_eigenvectors=False)
                    vals = np.sort(vals)[::-1]  # descending
                    return which, vals

                # which == 'LOW'  ---------------------------
                # Step 1: robust λ_min estimate
                lam_min_est = None
                try:
                    # largest algebraic of -S, then negate
                    vneg = eigsh(-S, k=1, which='LA', return_eigenvectors=False)[0]
                    lam_min_est = float(-vneg)
                except Exception:
                    try:
                        lam_min_est = float(eigsh(S, k=1, which='SA', return_eigenvectors=False)[0])
                    except Exception:
                        if n <= 1024:
                            lam_min_est = float(np.linalg.eigvalsh(S.toarray())[0])

                # Step 2: shift-invert around λ_min estimate to pull the actual k smallest
                if lam_min_est is not None:
                    try:
                        # a tiny negative epsilon so "nearest to shift" biases toward the minimum side
                        vals = eigsh(S, k=kk, sigma=lam_min_est - 1e-12, which='LM', return_eigenvectors=False)
                        vals = np.sort(vals)  # ascending
                        return which, vals
                    except Exception:
                        pass  # fall through to direct SA

                # Step 3: final fallback - ask for SA directly
                vals = eigsh(S, k=kk, which='SA', return_eigenvectors=False)
                vals = np.sort(vals)
                return which, vals

            else:
                # SVD path
                m, n = A.shape
                kk = int(min(k, max(1, min(m, n) - 1)))
                u, s, vt = svds(A.astype(np.float64), k=kk)
                s = np.sort(s)[::-1]
                return which, s

        worker = _Task(task)
        worker.finishedWithResult.connect(self._done)
        worker.error.connect(lambda e: self.log.appendPlainText(f"Error: {e}"))
        worker.start()
        self.worker = worker

    def _done(self, res: Tuple[str, np.ndarray]):
        which, arr = res
        self.table.setRowCount(0)
        for i, val in enumerate(arr):
            r = self.table.rowCount(); self.table.insertRow(r)
            self.table.setItem(r, 0, QtWidgets.QTableWidgetItem(str(i)))
            self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(f"{val:.8g}"))
        lbl = {'LOW':'Smallest eigenvalues (ascending)', 'HIGH':'Largest eigenvalues (descending)', 'SVD':'Top singular values (descending)'}[which]
        self.log.appendPlainText(f"Computed {lbl}: k={arr.size}")


class _Task(QtCore.QThread):
    finishedWithResult = QtCore.pyqtSignal(object)
    error = QtCore.pyqtSignal(str)
    def __init__(self, fn):
        super().__init__(); self.fn = fn
    def run(self):
        try:
            self.finishedWithResult.emit(self.fn())
        except Exception as e:
            self.error.emit(str(e))
