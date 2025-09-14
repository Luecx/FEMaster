
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