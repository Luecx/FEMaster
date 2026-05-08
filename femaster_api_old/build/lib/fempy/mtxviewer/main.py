# Project: Matrix Market Viewer (PyQt5 + pyqtgraph)
# Layout (files):
#   mtxviewer/
#     __init__.py
#     app.py            # entry point (run this)
#     matrix_view.py    # left matrix canvas + controls
#     tabs_overview.py  # Overview tab
#     tabs_spectral.py  # Spectral tab (table view)
#     tabs_decomp.py    # Decompositions & Solvers tab
#     tabs_formats.py   # Formats & memory tab
#     tabs_io.py        # I/O & export tab
#     style.py          # your dark stylesheet
#
# Save each block below into its filename. Then run:  python -m mtxviewer.app

# =============================
# File: mtxviewer/style.py
# =============================
STYLE = r"""
/* ---------- Base ---------- */
* { outline: 0; }
QWidget {
    background-color: #0f1115;
    color: #eef1f5;
    font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif; /* Segoe UI fallback */
    font-size: 15px;
}

/* Make toolbars, groupboxes, splitters consistent */
QToolBar { background: #0f1115; border-bottom: 1px solid #2b3340; }
QGroupBox { border: 1px solid #2b3340; border-radius: 10px; margin-top: 10px; background: #0f1115; }
QGroupBox::title { left: 10px; padding: 0 6px; color: #cbd5e1; background: transparent; }
QSplitter::handle { background: #141820; }

/* ---------- Buttons ---------- */
/* Qt doesn't support CSS classes; use a dynamic property instead. */
QPushButton[primary="true"] {
    background-color: #2563eb; color: white; border: none; border-radius: 10px; padding: 10px 14px; font-weight: 600;
}
QPushButton[primary="true"]:hover { background-color: #1d4ed8; }
QPushButton[primary="true"]:disabled { background-color: #1e293b; color: #8b9bb2; }
QPushButton { background: #141820; border: 1px solid #2b3340; border-radius: 8px; padding: 8px 12px; }
QPushButton:hover { background: #17202a; }
QToolButton { background: #141820; border: 1px solid #2b3340; border-radius: 8px; padding: 6px 10px; color: #cbd5e1; }
QToolButton:hover { background: #17202a; color: #fff; }

/* ---------- Inputs ---------- */
QComboBox, QSpinBox, QLineEdit { background: #141820; border: 1px solid #2b3340; border-radius: 8px; padding: 6px 10px; min-height: 30px; }
QComboBox QAbstractItemView { background: #0f131a; color: #eef1f5; border: 1px solid #2b3340; selection-background-color: #1f2937; }
QLabel { background: transparent; }

/* ---------- Tables / Views ---------- */
QTableWidget, QTreeView, QListView { background: #0f1315; border: 1px solid #2b3340; border-radius: 10px; gridline-color: #2b3340; selection-background-color: #1f2937; selection-color: #e5e7eb; }
QHeaderView::section { background: #121722; color: #e2e8f0; border: none; padding: 10px 12px; font-weight: 600; font-size: 14px; }
QTableCornerButton::section { background: #121722; border: none; }

/* ---------- Tabs ---------- */
QTabWidget::pane { border: 1px solid #2b3340; border-radius: 12px; padding: 8px; top: -2px; background: #0f1115; }
QTabBar::tab { background: #141820; color: #cbd5e1; padding: 10px 16px; margin-right: 6px; border: 1px solid #2b3340; border-bottom: 2px solid #2b3340; border-top-left-radius: 10px; border-top-right-radius: 10px; min-width: 140px; font-weight: 600; }
QTabBar::tab:hover { background: #17202a; color: #e5e7eb; }
QTabBar::tab:selected { background: #1a2230; color: #ffffff; border-bottom: 2px solid #2563eb; }

/* ---------- Scrollbars ---------- */
QScrollBar:vertical, QScrollBar:horizontal { background: #0f1115; border: none; margin: 0; }
QScrollBar:vertical { width: 10px; } QScrollBar:horizontal { height: 10px; }
QScrollBar::handle:vertical, QScrollBar::handle:horizontal { background: #2b3340; border-radius: 6px; min-height: 24px; min-width: 24px; }
QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover { background: #3a4453; }
QScrollBar::add-line, QScrollBar::sub-line, QScrollBar::add-page, QScrollBar::sub-page { background: transparent; border: none; }

/* ---------- Misc ---------- */
QPlainTextEdit, QTextEdit { background: #0f1315; border: 1px solid #2b3340; border-radius: 8px; }
"""

# =============================
# File: mtxviewer/matrix_view.py
# =============================
from typing import Optional, Tuple
import numpy as np
from scipy import sparse
from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg

class MatrixViewerWidget(QtWidgets.QWidget):
    modeChanged = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self._build_ui()

    def _build_ui(self):
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(0,0,0,0)
        # Local control bar (fixes your point C)
        ctrl = QtWidgets.QHBoxLayout()
        self.open_btn = QtWidgets.QPushButton("Open…")
        self.open_btn.setProperty('primary', True)
        self.view_combo = QtWidgets.QComboBox()
        self.view_combo.addItems(["scatter", "heatmap", "abs-heatmap", "signed-heatmap"])
        ctrl.addWidget(self.open_btn)
        ctrl.addSpacing(8)
        ctrl.addWidget(QtWidgets.QLabel("View:"))
        ctrl.addWidget(self.view_combo)
        ctrl.addStretch(1)
        v.addLayout(ctrl)

        # Graphics canvas
        self.gw = pg.GraphicsLayoutWidget()
        self.gw.setAntialiasing(False)
        self.view = self.gw.addViewBox()
        self.view.setAspectLocked(False)
        self.view.invertY(True)
        self.view.setMouseMode(self.view.RectMode)
        self.imgItem = pg.ImageItem(axisOrder='row-major')
        self.scatter = pg.ScatterPlotItem(size=2, brush=pg.mkBrush(230,230,230,180), pen=None)
        self.view.addItem(self.imgItem)
        self.view.addItem(self.scatter)
        v.addWidget(self.gw)

        self.mode = 'scatter'
        self.max_points = 1_000_000
        self.colormap = pg.colormap.get('inferno')
        self.imgItem.setLookupTable(self.colormap.getLookupTable(nPts=256))
        self.imgItem.setLevels((0.0, 1.0))  # tuple property, not callable

        self.view_combo.currentTextChanged.connect(self.set_mode)
        self.view.sigRangeChanged.connect(self._schedule_redraw)

    # -------- Public API --------
    def set_matrix(self, A: sparse.spmatrix):
        if not sparse.isspmatrix_csr(A):
            A = A.tocsr()
        self.A = A
        m, n = A.shape
        self.view.setRange(QtCore.QRectF(0, 0, n, m), padding=0)
        self._redraw(force=True)

    def set_mode(self, mode: str):
        if mode not in ('scatter','heatmap','abs-heatmap','signed-heatmap'):
            return
        self.mode = mode
        self.modeChanged.emit(mode)
        self._redraw(force=True)

    # -------- Drawing --------
    def _schedule_redraw(self):
        QtCore.QTimer.singleShot(40, lambda: self._redraw(False))

    def _redraw(self, force: bool):
        if self.A is None:
            return
        A = self.A
        m, n = A.shape
        rect = self.view.viewRect()
        x0 = max(int(np.floor(rect.left())), 0)
        y0 = max(int(np.floor(rect.top())), 0)
        x1 = min(int(np.ceil(rect.right())), n)
        y1 = min(int(np.ceil(rect.bottom())), m)
        if x1 <= x0 or y1 <= y0:
            return
        sub = A[y0:y1, x0:x1]

        if self.mode == 'scatter':
            self.imgItem.hide()
            rows, cols = sub.nonzero()
            pts = len(rows)
            if pts > self.max_points:
                idx = np.random.choice(pts, self.max_points, replace=False)
                rows, cols = rows[idx], cols[idx]
            spots = [{'pos': (x0 + int(c), y0 + int(r))} for r, c in zip(rows, cols)]
            self.scatter.setData(spots)
            self.scatter.show()
            return

        # Heatmaps
        self.scatter.hide()
        w = max(50, min(2000, int(rect.width())))
        h = max(50, min(2000, int(rect.height())))
        rows, cols = sub.nonzero()
        data = sub.data
        img = np.zeros((h, w), dtype=float)
        if rows.size:
            bx = np.clip(((cols - 0) / max(1, (x1 - x0)) * w).astype(int), 0, w-1)
            by = np.clip(((rows - 0) / max(1, (y1 - y0)) * h).astype(int), 0, h-1)
            if self.mode == 'heatmap':
                np.add.at(img, (by, bx), 1.0)
                levels = (float(img.min()), float(max(img.max(), 1.0)))
                self.imgItem.setLookupTable(self.colormap.getLookupTable(nPts=256))
                self.imgItem.setLevels(levels)
            elif self.mode == 'abs-heatmap':
                np.add.at(img, (by, bx), np.abs(data))
                levels = (float(img.min()), float(max(img.max(), 1e-12)))
                self.imgItem.setLookupTable(self.colormap.getLookupTable(nPts=256))
                self.imgItem.setLevels(levels)
            elif self.mode == 'signed-heatmap':
                np.add.at(img, (by, bx), data)
                mx = float(np.max(np.abs(img))) if img.size else 1.0
                self.imgItem.setLookupTable(pg.colormap.get('CET-D1').getLookupTable(nPts=256))
                self.imgItem.setLevels((-mx, mx))
        else:
            self.imgItem.setLookupTable(self.colormap.getLookupTable(nPts=256))
            self.imgItem.setLevels((0.0, 1.0))

        self.imgItem.setImage(img, levels=self.imgItem.levels)  # levels is a tuple property
        tr = QtGui.QTransform()
        sx = (x1 - x0) / max(1, w)
        sy = (y1 - y0) / max(1, h)
        tr.translate(x0, y0)
        tr.scale(sx, sy)
        self.imgItem.setTransform(tr)
        self.imgItem.show()

# =============================
# File: mtxviewer/tabs_overview.py
# =============================
from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
from scipy import sparse
from PyQt5 import QtWidgets

@dataclass
class MatrixInfo:
    shape: Tuple[int, int]
    nnz: int
    density: float
    is_square: bool
    is_symmetric_hint: bool
    fro_norm: float
    sym_abs_err: Optional[float]
    sym_rel_err: Optional[float]
    bandwidth: Optional[int]
    min_val: Optional[float]
    max_val: Optional[float]
    mean_val: Optional[float]
    std_val: Optional[float]
    eig_min: Optional[float]
    eig_max: Optional[float]
    spd: Optional[bool]

class OverviewTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._build()

    def _build(self):
        form = QtWidgets.QFormLayout(self)
        self.lbl_shape   = QtWidgets.QLabel("–")
        self.lbl_nnz     = QtWidgets.QLabel("–")
        self.lbl_density = QtWidgets.QLabel("–")
        self.lbl_sym     = QtWidgets.QLabel("–")
        self.lbl_sym_abs = QtWidgets.QLabel("–")
        self.lbl_sym_rel = QtWidgets.QLabel("–")
        self.lbl_bw      = QtWidgets.QLabel("–")
        self.lbl_min     = QtWidgets.QLabel("–")
        self.lbl_max     = QtWidgets.QLabel("–")
        self.lbl_mean    = QtWidgets.QLabel("–")
        self.lbl_std     = QtWidgets.QLabel("–")
        self.lbl_eigmin  = QtWidgets.QLabel("–")
        self.lbl_eigmax  = QtWidgets.QLabel("–")
        self.lbl_spd     = QtWidgets.QLabel("–")

        form.addRow("Shape", self.lbl_shape)
        form.addRow("NNZ", self.lbl_nnz)
        form.addRow("Density", self.lbl_density)
        form.addRow("Symmetric?", self.lbl_sym)
        form.addRow("Symmetry abs error ||A-A^T||_F", self.lbl_sym_abs)
        form.addRow("Symmetry rel error / ||A||_F", self.lbl_sym_rel)
        form.addRow("Bandwidth (max |i-j|)", self.lbl_bw)
        form.addRow("Min", self.lbl_min)
        form.addRow("Max", self.lbl_max)
        form.addRow("Mean", self.lbl_mean)
        form.addRow("Std", self.lbl_std)
        form.addRow("λ_min (sym part)", self.lbl_eigmin)
        form.addRow("λ_max (sym part)", self.lbl_eigmax)
        form.addRow("SPD?", self.lbl_spd)

    def set_info(self, info: MatrixInfo):
        self.lbl_shape.setText(f"{info.shape[0]} × {info.shape[1]}")
        self.lbl_nnz.setText(str(info.nnz))
        self.lbl_density.setText(f"{info.density*100:.6f}%")
        self.lbl_sym.setText("yes" if info.is_symmetric_hint else "no")
        self.lbl_sym_abs.setText("-" if info.sym_abs_err is None else f"{info.sym_abs_err:.6g}")
        self.lbl_sym_rel.setText("-" if info.sym_rel_err is None else f"{info.sym_rel_err:.6g}")
        self.lbl_bw.setText("-" if info.bandwidth is None else str(info.bandwidth))
        self.lbl_min.setText("-" if info.min_val is None else f"{info.min_val:.6g}")
        self.lbl_max.setText("-" if info.max_val is None else f"{info.max_val:.6g}")
        self.lbl_mean.setText("-" if info.mean_val is None else f"{info.mean_val:.6g}")
        self.lbl_std.setText("-" if info.std_val is None else f"{info.std_val:.6g}")
        self.lbl_eigmin.setText("-" if info.eig_min is None else f"{info.eig_min:.6g}")
        self.lbl_eigmax.setText("-" if info.eig_max is None else f"{info.eig_max:.6g}")
        if info.spd is None: self.lbl_spd.setText("-")
        else: self.lbl_spd.setText("yes" if info.spd else "no")

# =============================
# File: mtxviewer/tabs_formats.py
# =============================
import numpy as np
from scipy import sparse
from PyQt5 import QtWidgets

class FormatsTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Format", "dtype", "index dtype", "Bytes", "Human"])
        self.table.horizontalHeader().setStretchLastSection(True)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.table)

    def _mem_bytes(self, A: sparse.spmatrix, fmt: str, dtype: np.dtype, index_dtype=np.int32) -> int:
        m, n = A.shape
        nnz = A.nnz
        idx = np.dtype(index_dtype).itemsize
        val = np.dtype(dtype).itemsize
        if fmt in ('csr','csc'):
            ptr_len = m + 1 if fmt == 'csr' else n + 1
            return ptr_len*idx + nnz*(idx + val)
        elif fmt == 'coo':
            return nnz*(2*idx + val)
        return -1

    def _human(self, b: int) -> str:
        if b < 0: return "-"
        units = ['B','KB','MB','GB','TB']
        f = float(b); i=0
        while f >= 1024 and i < len(units)-1:
            f/=1024; i+=1
        return f"{f:.2f} {units[i]}"

    def set_matrix(self, A: sparse.spmatrix):
        self.table.setRowCount(0)
        dtypes = [np.float32, np.float64, np.complex64, np.complex128]
        idx_dtypes = [np.int32, np.int64]
        fmts = ['csr','csc','coo']
        for fmt in fmts:
            for dt in dtypes:
                for it in idx_dtypes:
                    r = self.table.rowCount(); self.table.insertRow(r)
                    bytes_ = self._mem_bytes(A, fmt, dt, it)
                    self.table.setItem(r,0, QtWidgets.QTableWidgetItem(fmt))
                    self.table.setItem(r,1, QtWidgets.QTableWidgetItem(np.dtype(dt).name))
                    self.table.setItem(r,2, QtWidgets.QTableWidgetItem(np.dtype(it).name))
                    self.table.setItem(r,3, QtWidgets.QTableWidgetItem(str(bytes_)))
                    self.table.setItem(r,4, QtWidgets.QTableWidgetItem(self._human(bytes_)))

# =============================
# File: mtxviewer/tabs_spectral.py
# =============================
from typing import Optional
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, svds
from PyQt5 import QtCore, QtWidgets

class SpectralTab(QtWidgets.QWidget):
    """
    Robust spectral routines for symmetric (use symmetric part S = (A+A^T)/2):
      - Smallest algebraic eigenvalues: eigsh(-S, which='LA') -> negate & ascend
      - Largest  algebraic eigenvalues: eigsh( S, which='LA') -> descend
      - Top singular values: svds(A)
    Results are presented in a table (index, value) with clear ordering.
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
        self.btn_low = QtWidgets.QPushButton("Smallest eigenvalues (algebraic)"); self.btn_low.setProperty('primary', True)
        self.btn_high= QtWidgets.QPushButton("Largest eigenvalues (algebraic)");  self.btn_high.setProperty('primary', True)
        self.btn_svd = QtWidgets.QPushButton("Top singular values");             self.btn_svd.setProperty('primary', True)
        ctl.addWidget(QtWidgets.QLabel("k:")); ctl.addWidget(self.k_spin); ctl.addStretch(1)
        ctl.addWidget(self.btn_low); ctl.addWidget(self.btn_high); ctl.addWidget(self.btn_svd)
        v.addLayout(ctl)
        self.table = QtWidgets.QTableWidget(0,2)
        self.table.setHorizontalHeaderLabels(["index","value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table)
        self.log = QtWidgets.QPlainTextEdit(); self.log.setReadOnly(True); self.log.setMaximumHeight(120)
        v.addWidget(self.log)
        self.btn_low.clicked.connect(lambda: self._compute('LOW'))
        self.btn_high.clicked.connect(lambda: self._compute('HIGH'))
        self.btn_svd.clicked.connect(lambda: self._compute('SVD'))

    def set_matrix(self, A: sparse.spmatrix):
        self.A = A
        self.table.setRowCount(0); self.log.setPlainText("")

    def _compute(self, which: str):
        if self.A is None: return
        k = int(self.k_spin.value())
        A = self.A
        def task():
            if which in ('LOW','HIGH'):
                # Symmetric part for eigenvalues
                S = ((A + A.T) * 0.5).asfptype()
                kk = min(k, max(1, min(S.shape)-2))
                if which == 'HIGH':
                    vals = eigsh(S, k=kk, which='LA', return_eigenvectors=False)
                    # Sort descending for clarity
                    vals = np.sort(vals)[::-1]
                    return which, vals
                else:  # LOW: smallest algebraic -> largest of -S, then negate
                    vals_neg = eigsh(-S, k=kk, which='LA', return_eigenvectors=False)
                    vals = -np.sort(vals_neg)[::-1]   # ascending smallest -> largest
                    return which, vals
            else:
                kk = min(k, max(1, min(A.shape)-1))
                # svds returns smallest to largest? We sort descending explicitly for "top"
                u, s, vt = svds(A, k=kk)
                s = np.sort(s)[::-1]
                return which, s
        worker = _Task(task)
        worker.finishedWithResult.connect(self._done)
        worker.error.connect(lambda e: self.log.appendPlainText(f"Error: {e}"))
        worker.start(); self.worker = worker

    def _done(self, res):
        which, arr = res
        self.table.setRowCount(0)
        for i, val in enumerate(arr):
            r = self.table.rowCount(); self.table.insertRow(r)
            self.table.setItem(r,0, QtWidgets.QTableWidgetItem(str(i)))
            self.table.setItem(r,1, QtWidgets.QTableWidgetItem(f"{val:.8g}"))
        label = { 'LOW': 'Smallest eigenvalues (ascending)',
                  'HIGH': 'Largest eigenvalues (descending)',
                  'SVD': 'Top singular values (descending)'}[which]
        self.log.appendPlainText(f"Computed {label}: k={arr.size}")

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

# =============================

from typing import Optional
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, svds
from PyQt5 import QtCore, QtWidgets

class SpectralTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self.worker: Optional[QtCore.QThread] = None
        self._build()

    def _build(self):
        v = QtWidgets.QVBoxLayout(self)
        ctl = QtWidgets.QHBoxLayout()
        self.k_spin = QtWidgets.QSpinBox(); self.k_spin.setRange(1, 500); self.k_spin.setValue(10)
        self.btn_low = QtWidgets.QPushButton("Smallest eigenvalues (SA)"); self.btn_low.setProperty('primary', True)
        self.btn_high= QtWidgets.QPushButton("Largest eigenvalues (LA)"); self.btn_high.setProperty('primary', True)
        self.btn_svd = QtWidgets.QPushButton("Top singular values"); self.btn_svd.setProperty('primary', True)
        ctl.addWidget(QtWidgets.QLabel("k:")); ctl.addWidget(self.k_spin); ctl.addStretch(1)
        ctl.addWidget(self.btn_low); ctl.addWidget(self.btn_high); ctl.addWidget(self.btn_svd)
        v.addLayout(ctl)
        self.table = QtWidgets.QTableWidget(0,2)
        self.table.setHorizontalHeaderLabels(["index","value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table)
        self.log = QtWidgets.QPlainTextEdit(); self.log.setReadOnly(True); self.log.setMaximumHeight(120)
        v.addWidget(self.log)
        self.btn_low.clicked.connect(lambda: self._compute('SA'))
        self.btn_high.clicked.connect(lambda: self._compute('LA'))
        self.btn_svd.clicked.connect(lambda: self._compute('SVD'))

    def set_matrix(self, A: sparse.spmatrix):
        self.A = A
        self.table.setRowCount(0); self.log.setPlainText("")

    def _compute(self, which: str):
        if self.A is None: return
        k = int(self.k_spin.value())
        A = self.A
        def task():
            if which in ('SA','LA'):
                # Use symmetric part to avoid surprises
                S = (A + A.T) * 0.5
                kk = min(k, min(S.shape)-2)
                vals = eigsh(S, k=kk, which=which, return_eigenvectors=False)
                vals.sort()
                return which, vals
            else:
                kk = min(k, min(A.shape)-1)
                u, s, vt = svds(A, k=kk)
                s.sort()
                return which, s
        worker = _Task(task)
        worker.finishedWithResult.connect(self._done)
        worker.error.connect(lambda e: self.log.appendPlainText(f"Error: {e}"))
        worker.start(); self.worker = worker

    def _done(self, res):
        which, arr = res
        self.table.setRowCount(0)
        for i, val in enumerate(arr):
            r = self.table.rowCount(); self.table.insertRow(r)
            self.table.setItem(r,0, QtWidgets.QTableWidgetItem(str(i)))
            self.table.setItem(r,1, QtWidgets.QTableWidgetItem(f"{val:.8g}"))
        self.log.appendPlainText(f"Computed {which}: k={arr.size}")

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

# =============================
# File: mtxviewer/tabs_decomp.py
# =============================
from typing import Optional
import numpy as np
from scipy import sparse
from PyQt5 import QtWidgets

class DecompTab(QtWidgets.QWidget):
    """Show which solvers/decompositions are applicable (green/red)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self.spd: Optional[bool] = None
        self.is_square: bool = False
        self._build()

    def _build(self):
        v = QtWidgets.QVBoxLayout(self)
        self.table = QtWidgets.QTableWidget(0,3)
        self.table.setHorizontalHeaderLabels(["Method", "Type", "Applicable?"])
        self.table.horizontalHeader().setStretchLastSection(True)
        v.addWidget(self.table)
        self.note = QtWidgets.QLabel("Applicability is inferred from square/symmetry/SPD checks.")
        self.note.setWordWrap(True)
        v.addWidget(self.note)

    def set_matrix(self, A: sparse.spmatrix, spd: Optional[bool]):
        self.A = A; self.is_square = (A.shape[0]==A.shape[1]); self.spd = spd
        self._refresh()

    def _add_row(self, name, typ, ok):
        r = self.table.rowCount(); self.table.insertRow(r)
        self.table.setItem(r,0, QtWidgets.QTableWidgetItem(name))
        self.table.setItem(r,1, QtWidgets.QTableWidgetItem(typ))
        cell = QtWidgets.QTableWidgetItem("✔" if ok else "✖")
        cell.setForeground(QtWidgets.QBrush(QtWidgets.QColor("#86efac" if ok else "#fca5a5")))
        self.table.setItem(r,2, cell)

    def _refresh(self):
        self.table.setRowCount(0)
        sq = self.is_square
        spd = bool(self.spd) if self.spd is not None else False
        # Direct
        self._add_row("Cholesky", "Direct", sq and spd)
        self._add_row("LDL^T (symmetric indefinite)", "Direct", sq and self.spd is not None)
        self._add_row("LU (general)", "Direct", sq)
        # Iterative
        self._add_row("CG (SPD)", "Iterative", sq and spd)
        self._add_row("MINRES (symmetric)", "Iterative", sq and self.spd is not None)
        self._add_row("GMRES (general)", "Iterative", sq)
        self._add_row("BiCGSTAB (general)", "Iterative", sq)

# =============================
# File: mtxviewer/tabs_io.py
# =============================
from typing import Optional
from PyQt5 import QtWidgets
from scipy import sparse
import pyqtgraph as pg

class IOExportTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        v = QtWidgets.QVBoxLayout(self)
        self.btn_npz = QtWidgets.QPushButton("Save as .npz (CSR)"); self.btn_npz.setProperty('primary', True)
        self.btn_png = QtWidgets.QPushButton("Save current view (PNG)"); self.btn_png.setProperty('primary', True)
        v.addWidget(self.btn_npz); v.addWidget(self.btn_png)
        self.out = QtWidgets.QPlainTextEdit(); self.out.setReadOnly(True); v.addWidget(self.out)
        self.viewer = None; self.A: Optional[sparse.spmatrix] = None
        self.btn_npz.clicked.connect(self._save_npz)
        self.btn_png.clicked.connect(self._save_png)

    def link(self, viewer, A: Optional[sparse.spmatrix]):
        self.viewer = viewer; self.A = A

    def _save_npz(self):
        if self.A is None: return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save .npz", filter="NumPy zip (*.npz)")
        if not path: return
        try:
            sparse.save_npz(path, self.A.tocsr()); self.out.appendPlainText(f"Saved to {path}")
        except Exception as e:
            self.out.appendPlainText(f"Error: {e}")

    def _save_png(self):
        if self.viewer is None: return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save PNG", filter="PNG Image (*.png)")
        if not path: return
        try:
            exporter = pg.exporters.ImageExporter(self.viewer.gw.scene())
            exporter.parameters()['width'] = 2000
            exporter.export(path)
            self.out.appendPlainText(f"Snapshot saved to {path}")
        except Exception as e:
            self.out.appendPlainText(f"Error: {e}")

# =============================
# File: mtxviewer/app.py
# =============================
import os
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
        super().__init__(); self.path = path
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
        M = mmread(path)
        return M.tocsr() if sparse.issparse(M) else sparse.csr_matrix(M)

class EigBoundsTask(QtCore.QThread):
    finishedWithResult = QtCore.pyqtSignal(object)
    error = QtCore.pyqtSignal(str)
    def __init__(self, A: sparse.spmatrix):
        super().__init__(); self.A = A
    def run(self):
        try:
            if self.A.shape[0] != self.A.shape[1]:
                self.finishedWithResult.emit((None, None, None))
                return
            S = (self.A + self.A.T) * 0.5
            # Robust: λ_max(A) = LA(S), λ_min(A) = -LA(-S)
            lmax = float(eigsh(S, k=1, which='LA', return_eigenvectors=False)[0])
            lmin = float(-eigsh(-S, k=1, which='LA', return_eigenvectors=False)[0])
            spd = bool(lmin > 0)
            self.finishedWithResult.emit((lmin, lmax, spd))
        except Exception as e:
            self.error.emit(str(e))

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Matrix Market Viewer — Sparse & Spectral")
        self.resize(1480, 960)
        pg.setConfigOption('background', '#0f1115'); pg.setConfigOption('foreground', 'w')

        splitter = QtWidgets.QSplitter(); splitter.setOrientation(QtCore.Qt.Horizontal)
        # Left panel (controls + viewer)
        self.viewer = MatrixViewerWidget()
        splitter.addWidget(self.viewer)
        # Right tabs
        self.tabs = QtWidgets.QTabWidget()
        self.tab_overview = OverviewTab()
        self.tab_spectral = SpectralTab()
        self.tab_decomp   = DecompTab()
        self.tab_formats  = FormatsTab()
        self.tab_io       = IOExportTab()
        for name, w in [("Overview", self.tab_overview), ("Spectral", self.tab_spectral), ("Decompositions", self.tab_decomp), ("Formats & Memory", self.tab_formats), ("I/O & Export", self.tab_io)]:
            self.tabs.addTab(w, name)
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(0, 3); splitter.setStretchFactor(1, 2)
        self.setCentralWidget(splitter)

        # Wire open
        self.viewer.open_btn.clicked.connect(self._open)

        self.A: Optional[sparse.spmatrix] = None

    # ---- helpers ----
    def _build_info(self, A: sparse.spmatrix, eigs: Tuple[Optional[float],Optional[float],Optional[bool]]=(None,None,None)) -> MatrixInfo:
        m, n = A.shape
        is_sym_hint = (m==n) and (A.nnz == A.T.nnz)
        if A.nnz:
            mn = float(A.data.min()); mx = float(A.data.max()); mean=float(A.data.mean()); std=float(A.data.std())
        else:
            mn=mx=mean=std=None
        fro = float(np.sqrt((A.power(2)).sum())) if A.nnz else 0.0
        sym_abs = sym_rel = None
        try:
            D = A - A.T
            sym_abs = float(np.sqrt(D.power(2).sum()))
            sym_rel = float(sym_abs / max(fro, 1e-30))
        except Exception:
            pass
        try:
            r,c = A.nonzero(); bw = int(np.max(np.abs(r-c))) if r.size else 0
        except Exception:
            bw = None
        lmin, lmax, spd = eigs
        return MatrixInfo(
            shape=(m,n), nnz=A.nnz, density=A.nnz/max(1,m*n), is_square=(m==n), is_symmetric_hint=is_sym_hint,
            fro_norm=fro, sym_abs_err=sym_abs, sym_rel_err=sym_rel, bandwidth=bw,
            min_val=mn, max_val=mx, mean_val=mean, std_val=std, eig_min=lmin, eig_max=lmax, spd=spd
        )

    def _open(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open matrix", filter="Matrix files (*.mtx *.mtx.gz *.npz *.npy);;All files (*)")
        if not path: return
        loader = LoaderTask(path)
        loader.finishedWithResult.connect(self._loaded)
        loader.error.connect(lambda e: QtWidgets.QMessageBox.critical(self, "Load error", e))
        loader.start(); self._loader = loader

    def _loaded(self, A: sparse.spmatrix):
        self.A = A.tocsr()
        self.viewer.set_matrix(self.A)
        if self.A.shape[0]==self.A.shape[1]:
            eigtask = EigBoundsTask(self.A)
            eigtask.finishedWithResult.connect(self._post_eigs)
            eigtask.error.connect(lambda e: self._post_eigs((None,None,None)))
            eigtask.start(); self._eig = eigtask
        else:
            info = self._build_info(self.A, (None,None,None))
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
    app = QtWidgets.QApplication([])
    app.setStyleSheet(STYLE)
    win = MainWindow(); win.show()
    return app.exec_()

if __name__ == '__main__':
    main()


# =============================
# File: mtxviewer/__init__.py
# =============================
# empty init to make it a package
