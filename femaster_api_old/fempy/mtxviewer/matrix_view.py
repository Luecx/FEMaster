# =============================
# File: mtxviewer/matrix_view.py
# =============================
from typing import Optional, Tuple
import numpy as np
from scipy import sparse
from PyQt5 import QtCore, QtWidgets, QtGui
import pyqtgraph as pg


class _GridOverlay(pg.GraphicsObject):
    """
    Cell-aligned grid overlay with user-configurable major/minor logic.

    Zoom levels (from pixel size per cell):
      level 0 (finest):   minor = 1,                   major = base
      level L >= 1:       minor = base*growth^(L-1),   major = base*growth^L

    Always shows both minor & major for the current level. Lines are anchored to 0,
    so they do not drift while panning.
    """
    def __init__(self, view: pg.ViewBox, ncols: int = 0, nrows: int = 0):
        super().__init__()
        self.view = view
        self.ncols = int(ncols)
        self.nrows = int(nrows)
        self.enabled = True
        self.base_major = 5
        self.growth_factor = 2
        self.min_px_per_cell = 2.0  # below this we don't draw any grid
        self.setZValue(1000)

    def set_enabled(self, on: bool):
        self.enabled = bool(on)
        self.update()

    def set_extent(self, ncols: int, nrows: int):
        self.prepareGeometryChange()
        self.ncols = int(ncols)
        self.nrows = int(nrows)
        self.update()

    def set_base_major(self, base: int):
        self.base_major = max(1, int(base))
        self.update()

    def set_growth_factor(self, g: int):
        self.growth_factor = max(1, int(g))
        self.update()

    def boundingRect(self) -> QtCore.QRectF:
        return QtCore.QRectF(0.0, 0.0, float(self.ncols), float(self.nrows))

    def _zoom_level(self, px_per_cell: float) -> Optional[int]:
        if px_per_cell < self.min_px_per_cell:
            return None
        if px_per_cell >= 16: return 0
        if px_per_cell >= 8:  return 1
        if px_per_cell >= 4:  return 2
        return 3

    def _strides_for_level(self, level: int) -> Tuple[int, int]:
        base, g = self.base_major, self.growth_factor
        if level <= 0:
            return 1, max(1, base)
        minor = max(1, int(base * (g ** (level - 1))))
        major = max(minor + 1, int(base * (g ** level)))
        if major % minor != 0:  # clean alignment
            major = ((major // minor) + 1) * minor
        return minor, major

    def paint(self, painter: QtGui.QPainter, opt, widget=None):
        if not self.enabled or self.ncols <= 0 or self.nrows <= 0:
            return

        psize = self.view.viewPixelSize()
        if not psize or psize[0] <= 0 or psize[1] <= 0:
            return
        px_per_cell_x = 1.0 / float(psize[0])
        px_per_cell_y = 1.0 / float(psize[1])

        level_x = self._zoom_level(px_per_cell_x)
        level_y = self._zoom_level(px_per_cell_y)
        if level_x is None or level_y is None:
            return

        minor_x, major_x = self._strides_for_level(level_x)
        minor_y, major_y = self._strides_for_level(level_y)

        # Visible rect (clipped)
        vr = self.view.viewRect()
        x0 = max(int(np.floor(max(0.0, vr.left()))), 0)
        x1 = min(int(np.ceil(min(float(self.ncols), vr.right()))), self.ncols)
        y0 = max(int(np.floor(max(0.0, vr.top()))), 0)
        y1 = min(int(np.ceil(min(float(self.nrows), vr.bottom()))), self.nrows)
        if x1 <= x0 or y1 <= y0:
            return

        MAX_LINES = 4000  # total budget

        def positions(a0, a1, stride):
            if stride < 1: stride = 1
            first = ((a0 + stride - 1) // stride) * stride  # anchor to 0 multiples
            return np.arange(first, a1 + 1, stride, dtype=int)

        # Throttle while keeping anchoring (increase stride by integer factor)
        def throttle(a0, a1, stride, budget):
            s = max(1, int(stride))
            xs = positions(a0, a1, s)
            while xs.size > budget:
                s *= 2
                xs = positions(a0, a1, s)
            return s, xs

        # Ensure major is multiple of minor before throttling
        if major_x % minor_x != 0:
            major_x = ((major_x // minor_x) + 1) * minor_x
        if major_y % minor_y != 0:
            major_y = ((major_y // minor_y) + 1) * minor_y

        # Throttle minors & majors independently but preserve multiples
        minor_x, xs_minor = throttle(x0, x1, minor_x, MAX_LINES // 2)
        if major_x % minor_x != 0:
            major_x = ((major_x // minor_x) + 1) * minor_x
        major_x, xs_major = throttle(x0, x1, major_x, MAX_LINES // 2)

        minor_y, ys_minor = throttle(y0, y1, minor_y, MAX_LINES // 2)
        if major_y % minor_y != 0:
            major_y = ((major_y // minor_y) + 1) * minor_y
        major_y, ys_major = throttle(y0, y1, major_y, MAX_LINES // 2)

        # Remove majors from minors to avoid overdrawing
        if xs_minor.size and xs_major.size:
            xs_minor = xs_minor[~np.isin(xs_minor, xs_major)]
        if ys_minor.size and ys_major.size:
            ys_minor = ys_minor[~np.isin(ys_minor, ys_major)]

        if xs_minor.size == 0 and xs_major.size == 0 and ys_minor.size == 0 and ys_major.size == 0:
            return

        # Clip endpoints to current view
        ylo, yhi = (float(min(vr.top(), vr.bottom())), float(max(vr.top(), vr.bottom())))
        xlo, xhi = (float(min(vr.left(), vr.right())), float(max(vr.left(), vr.right())))
        ylo = max(0.0, ylo); yhi = min(float(self.nrows), yhi)
        xlo = max(0.0, xlo); xhi = min(float(self.ncols), xhi)

        # Pens
        pal = widget.palette() if widget else QtWidgets.QApplication.palette()
        base_col = pal.color(QtGui.QPalette.WindowText)
        minor_color = QtGui.QColor(base_col); minor_color.setAlpha(60)
        major_color = QtGui.QColor(base_col); major_color.setAlpha(130)
        pen_minor = QtGui.QPen(minor_color); pen_minor.setWidth(0)
        pen_major = QtGui.QPen(major_color); pen_major.setWidth(0)

        painter.setRenderHint(QtGui.QPainter.Antialiasing, False)

        # Draw verticals
        if xs_minor.size:
            painter.setPen(pen_minor)
            painter.drawLines([QtCore.QLineF(float(x), ylo, float(x), yhi) for x in xs_minor])
        if xs_major.size:
            painter.setPen(pen_major)
            painter.drawLines([QtCore.QLineF(float(x), ylo, float(x), yhi) for x in xs_major])

        # Draw horizontals
        if ys_minor.size:
            painter.setPen(pen_minor)
            painter.drawLines([QtCore.QLineF(xlo, float(y), xhi, float(y)) for y in ys_minor])
        if ys_major.size:
            painter.setPen(pen_major)
            painter.drawLines([QtCore.QLineF(xlo, float(y), xhi, float(y)) for y in ys_major])


class MatrixViewerWidget(QtWidgets.QWidget):
    modeChanged = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self.current_img: Optional[np.ndarray] = None
        self._build_ui()
        self.MAX_TOTAL_PIXELS = 16_000_000
        self.MIN_SIDE = 512

    def _build_ui(self):
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)

        ctrl = QtWidgets.QHBoxLayout()
        self.open_btn = QtWidgets.QPushButton("Open…")
        self.open_btn.setProperty("primary", True)

        self.view_combo = QtWidgets.QComboBox()
        self.view_combo.addItems(["count heatmap", "abs-heatmap", "signed-heatmap"])

        self.grid_cb = QtWidgets.QCheckBox("Grid"); self.grid_cb.setChecked(True)
        self.grid_base = QtWidgets.QSpinBox(); self.grid_base.setRange(1, 999); self.grid_base.setValue(5)
        self.grid_base.setToolTip("Major grid spacing at finest zoom (cells)")
        self.grid_growth = QtWidgets.QSpinBox(); self.grid_growth.setRange(1, 16); self.grid_growth.setValue(2)
        self.grid_growth.setToolTip("Multiply major spacing as you zoom out (e.g., 4 → 3, 12, 48, …)")

        ctrl.addWidget(self.open_btn)
        ctrl.addSpacing(8)
        ctrl.addWidget(QtWidgets.QLabel("View:"))
        ctrl.addWidget(self.view_combo)
        ctrl.addSpacing(12)
        ctrl.addWidget(self.grid_cb)
        ctrl.addWidget(QtWidgets.QLabel("Base:"))
        ctrl.addWidget(self.grid_base)
        ctrl.addWidget(QtWidgets.QLabel("Growth ×"))
        ctrl.addWidget(self.grid_growth)
        ctrl.addStretch(1)
        v.addLayout(ctrl)

        self.gw = pg.GraphicsLayoutWidget()
        self.gw.setAntialiasing(False)
        self.view = self.gw.addViewBox(lockAspect=True, enableMenu=False)
        self.view.invertY(True)

        win_color = self.palette().color(QtGui.QPalette.Window)
        self.gw.setBackground(pg.mkColor(win_color))
        if hasattr(self.view, "setBackgroundColor"):
            self.view.setBackgroundColor(pg.mkColor(win_color))

        self.imgItem = pg.ImageItem(axisOrder="row-major")
        self.imgItem.setAutoDownsample(True)
        self.view.addItem(self.imgItem)

        self.grid = _GridOverlay(self.view)
        self.view.addItem(self.grid)

        v.addWidget(self.gw)

        self.cmap_seq = pg.colormap.get("inferno")
        self.cmap_div = pg.colormap.get("CET-D1")

        self.view_combo.currentTextChanged.connect(self.set_mode)
        self.grid_cb.toggled.connect(self._toggle_grid)
        self.grid_base.valueChanged.connect(self._update_grid_params)
        self.grid_growth.valueChanged.connect(self._update_grid_params)
        self.view.sigRangeChanged.connect(lambda *_: self.grid.update())

    def set_matrix(self, A: sparse.spmatrix):
        if not sparse.isspmatrix_csr(A):
            A = A.tocsr()
        self.A = A
        self.set_mode(self.view_combo.currentText())

    def set_mode(self, mode: str):
        if self.A is None:
            return
        self.mode = mode
        self.modeChanged.emit(mode)
        self._render_global_image()

    def _toggle_grid(self, on: bool):
        self.grid.set_enabled(on)
        self.grid.update()

    def _update_grid_params(self):
        self.grid.set_base_major(self.grid_base.value())
        self.grid.set_growth_factor(self.grid_growth.value())
        self.grid.update()

    def _target_size(self, shape: Tuple[int, int]) -> Tuple[int, int]:
        m, n = shape
        total = int(m) * int(n)
        if total <= self.MAX_TOTAL_PIXELS:
            return int(n), int(m)
        aspect = float(n) / max(1.0, float(m))
        w = int(np.sqrt(self.MAX_TOTAL_PIXELS * aspect))
        h = max(1, int(w / max(aspect, 1e-12)))
        w = max(self.MIN_SIDE, w)
        h = max(self.MIN_SIDE, h)
        return int(w), int(h)

    def _render_global_image(self):
        A = self.A
        m, n = A.shape
        coo = A.tocoo()
        rows, cols, data = coo.row, coo.col, coo.data

        out_w, out_h = self._target_size((m, n))
        img = np.zeros((out_h, out_w), dtype=np.float32)

        if rows.size:
            bx = (cols.astype(np.float64) * (out_w / max(1.0, float(n)))).astype(np.int64)
            by = (rows.astype(np.float64) * (out_h / max(1.0, float(m)))).astype(np.int64)
            np.clip(bx, 0, out_w - 1, out=bx)
            np.clip(by, 0, out_h - 1, out=by)

            if self.mode == "count heatmap":
                np.add.at(img, (by, bx), 1.0)
                lut = self.cmap_seq.getLookupTable(nPts=256)
                vmax = float(img.max()) if img.size else 1.0
                self.imgItem.setLookupTable(lut)
                self.imgItem.setLevels((0.0, max(vmax, 1.0)))
            elif self.mode == "abs-heatmap":
                np.add.at(img, (by, bx), np.abs(data).astype(np.float32))
                lut = self.cmap_seq.getLookupTable(nPts=256)
                vmax = float(img.max()) if img.size else 1.0
                self.imgItem.setLookupTable(lut)
                self.imgItem.setLevels((0.0, max(vmax, 1e-12)))
            else:
                np.add.at(img, (by, bx), data.astype(np.float32))
                mx = float(np.max(np.abs(img))) if img.size else 1.0
                lut = self.cmap_div.getLookupTable(nPts=256)
                self.imgItem.setLookupTable(lut)
                self.imgItem.setLevels((-mx, mx))
        else:
            lut = self.cmap_seq.getLookupTable(nPts=256)
            self.imgItem.setLookupTable(lut)
            self.imgItem.setLevels((0.0, 1.0))

        self.current_img = img
        self.imgItem.setImage(img, autoLevels=False)
        self.imgItem.setRect(QtCore.QRectF(0, 0, float(n), float(m)))
        self.view.setRange(QtCore.QRectF(0, 0, float(n), float(m)), padding=0)

        self.grid.set_extent(n, m)
        self._update_grid_params()
        self.grid.update()
