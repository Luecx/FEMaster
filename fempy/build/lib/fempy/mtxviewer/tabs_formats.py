
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
