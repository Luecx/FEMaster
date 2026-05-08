# =============================
# File: mtxviewer/app.py
# =============================
import os
import sys
import argparse
from typing import Optional, Tuple
import numpy as np
from scipy import sparse
from scipy.io import mmread
from scipy.sparse.linalg import eigsh
from PyQt5 import QtCore, QtWidgets
import pyqtgraph as pg

from .style import STYLE
from .matrix_view import MatrixViewerWidget
from .tabs_overview import OverviewTab, MatrixInfo
from .tabs_spectral import SpectralTab
from .tabs_decomp import DecompTab
from .tabs_formats import FormatsTab
from .tabs_io import IOExportTab


class LoaderTask(QtCore.QThread):
    finishedWithResult = QtCore.pyqtSignal(object)
    error = QtCore.pyqtSignal(str)

    def __init__(self, path: str):
        super().__init__()
        self.path = path

    def run(self):
        try:
            A = self._load_matrix(self.path)
            self.finishedWithResult.emit(A)
        except Exception as e:
            self.error.emit(str(e))

    def _load_matrix(self, path: str) -> sparse.spmatrix:
        p = path.lower()
        if p.endswith('.mtx') or p.endswith('.mtx.gz'):
            M = mmread(path)
            return M.tocsr() if sparse.issparse(M) else sparse.csr_matrix(M)
        if p.endswith('.npz'):
            return sparse.load_npz(path).tocsr()
        if p.endswith('.npy'):
            return sparse.csr_matrix(np.load(path))
        # Fallback: try MatrixMarket
        M = mmread(path)
        return M.tocsr() if sparse.issparse(M) else sparse.csr_matrix(M)


class EigBoundsTask(QtCore.QThread):
    finishedWithResult = QtCore.pyqtSignal(object)
    error = QtCore.pyqtSignal(str)

    def __init__(self, A: sparse.spmatrix):
        super().__init__()
        self.A = A

    def run(self):
        try:
            m, n = self.A.shape
            if m != n:
                self.finishedWithResult.emit((None, None, None))
                return

            # Work on symmetric part, float64 CSC for ARPACK stability
            S = ((self.A + self.A.T) * 0.5).astype(np.float64).tocsc()
            n = S.shape[0]

            # Helper: safe largest-eig (LA) with small fallbacks
            def largest_eig(X):
                try:
                    val = float(eigsh(X, k=1, which='LA', return_eigenvectors=False)[0])
                    return val
                except Exception:
                    # Tiny dense fallback
                    if n <= 1024:
                        return float(np.linalg.eigvalsh(X.toarray())[-1])
                    raise

            # λ_max(S)
            lmax = largest_eig(S)

            # λ_min(S) = - λ_max(-S)
            try:
                lmin = -largest_eig(-S)
            except Exception:
                # Fallback: try SA directly
                try:
                    lmin = float(eigsh(S, k=1, which='SA', return_eigenvectors=False)[0])
                except Exception:
                    # Last resort: dense for very small n, else give up gracefully
                    if n <= 1024:
                        lmin = float(np.linalg.eigvalsh(S.toarray())[0])
                    else:
                        self.finishedWithResult.emit((None, lmax, None))
                        return

            spd = bool(lmin > 0.0)
            self.finishedWithResult.emit((lmin, lmax, spd))

        except Exception as e:
            self.error.emit(str(e))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, initial_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("Matrix Market Viewer — Sparse & Spectral")
        self.resize(1480, 960)
        pg.setConfigOption('background', '#0f1115')
        pg.setConfigOption('foreground', 'w')

        splitter = QtWidgets.QSplitter()
        splitter.setOrientation(QtCore.Qt.Horizontal)

        # Left panel (controls + viewer)
        self.viewer = MatrixViewerWidget()
        splitter.addWidget(self.viewer)

        # Right tabs
        self.tabs = QtWidgets.QTabWidget()
        self.tab_overview = OverviewTab()
        self.tab_spectral = SpectralTab()
        self.tab_decomp = DecompTab()
        self.tab_formats = FormatsTab()
        self.tab_io = IOExportTab()
        for name, w in [
            ("Overview", self.tab_overview),
            ("Spectral", self.tab_spectral),
            ("Decompositions", self.tab_decomp),
            ("Formats & Memory", self.tab_formats),
            ("I/O & Export", self.tab_io),
        ]:
            self.tabs.addTab(w, name)

        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        self.setCentralWidget(splitter)

        # Wire open
        self.viewer.open_btn.clicked.connect(self._open)

        self.A: Optional[sparse.spmatrix] = None

        # Auto-load if a path was provided
        if initial_path:
            self._open_path(initial_path)

    # ---- helpers ----
    def _build_info(
            self,
            A: sparse.spmatrix,
            eigs: Tuple[Optional[float], Optional[float], Optional[bool]] = (None, None, None),
    ) -> MatrixInfo:
        m, n = A.shape
        is_sym_hint = (m == n) and (A.nnz == A.T.nnz)
        if A.nnz:
            mn = float(A.data.min())
            mx = float(A.data.max())
            mean = float(A.data.mean())
            std = float(A.data.std())
        else:
            mn = mx = mean = std = None
        fro = float(np.sqrt((A.power(2)).sum())) if A.nnz else 0.0
        sym_abs = sym_rel = None
        try:
            D = A - A.T
            sym_abs = float(np.sqrt(D.power(2).sum()))
            sym_rel = float(sym_abs / max(fro, 1e-30))
        except Exception:
            pass
        try:
            r, c = A.nonzero()
            bw = int(np.max(np.abs(r - c))) if r.size else 0
        except Exception:
            bw = None
        lmin, lmax, spd = eigs
        return MatrixInfo(
            shape=(m, n),
            nnz=A.nnz,
            density=A.nnz / max(1, m * n),
            is_square=(m == n),
            is_symmetric_hint=is_sym_hint,
            fro_norm=fro,
            sym_abs_err=sym_abs,
            sym_rel_err=sym_rel,
            bandwidth=bw,
            min_val=mn,
            max_val=mx,
            mean_val=mean,
            std_val=std,
            eig_min=lmin,
            eig_max=lmax,
            spd=spd,
        )

    def _open(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open matrix",
            filter="Matrix files (*.mtx *.mtx.gz *.npz *.npy);;All files (*)",
        )
        if not path:
            return
        self._open_path(path)

    def _open_path(self, path: str):
        path = os.path.abspath(os.path.expanduser(path))
        if not os.path.exists(path):
            QtWidgets.QMessageBox.critical(self, "Load error", f"File not found:\n{path}")
            return
        loader = LoaderTask(path)
        loader.finishedWithResult.connect(self._loaded)
        loader.error.connect(lambda e: QtWidgets.QMessageBox.critical(self, "Load error", e))
        loader.start()
        self._loader = loader

    def _loaded(self, A: sparse.spmatrix):
        self.A = A.tocsr()
        self.viewer.set_matrix(self.A)
        if self.A.shape[0] == self.A.shape[1]:
            eigtask = EigBoundsTask(self.A)
            eigtask.finishedWithResult.connect(self._post_eigs)
            eigtask.error.connect(lambda e: self._post_eigs((None, None, None)))
            eigtask.start()
            self._eig = eigtask
        else:
            info = self._build_info(self.A, (None, None, None))
            self._push_info(info)

    def _post_eigs(self, triple):
        info = self._build_info(self.A, triple)
        self._push_info(info)

    def _push_info(self, info: MatrixInfo):
        self.tab_overview.set_info(info)
        self.tab_formats.set_matrix(self.A)
        self.tab_spectral.set_matrix(self.A)
        self.tab_decomp.set_matrix(self.A, info.spd)
        self.tab_io.link(self.viewer, self.A)


# ---- app entry ----
def main():
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

    # Parse optional file path
    parser = argparse.ArgumentParser(description="Matrix Market Viewer")
    parser.add_argument("path", nargs="?", help="Matrix file to open (.mtx/.mtx.gz/.npz/.npy)")
    parser.add_argument("--file", "-f", dest="path_flag", help="Matrix file to open")
    args = parser.parse_args()

    initial_path = args.path_flag or args.path

    app = QtWidgets.QApplication(sys.argv[:1])  # let argparse consume extras
    app.setStyleSheet(STYLE)
    win = MainWindow(initial_path=initial_path)
    win.show()
    return app.exec_()


if __name__ == "__main__":
    main()
