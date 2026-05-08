
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