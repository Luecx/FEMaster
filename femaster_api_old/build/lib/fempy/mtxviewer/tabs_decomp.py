# =============================
# File: mtxviewer/tabs_decomp.py
# =============================
from typing import Optional
import numpy as np
from scipy import sparse
from PyQt5 import QtWidgets, QtGui, QtCore


class DecompTab(QtWidgets.QWidget):
    """
    Show which solvers/decompositions are applicable for the loaded matrix.

    Heuristics used:
      - square        := (m == n)
      - sym_hint      := square and (A.nnz == A.T.nnz)   (cheap structural symmetry hint)
      - spd           := from Overview (λ_min > 0). May be None if unknown.
      - rectangular   := not square

    Columns:
      Method | Type | Requirements | Applicable?
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.A: Optional[sparse.spmatrix] = None
        self.spd: Optional[bool] = None
        self._build()

    # ---------- UI ----------
    def _build(self):
        v = QtWidgets.QVBoxLayout(self)

        self.table = QtWidgets.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Method", "Type", "Requirements", "Applicable?"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
        v.addWidget(self.table)

        self.note = QtWidgets.QLabel(
            "Applicability is inferred from shape/symmetry/SPD checks.\n"
            "• sym_hint is a cheap structural proxy. Exact symmetry is verified elsewhere.\n"
            "• Least-squares solvers apply when the system is over/underdetermined."
        )
        self.note.setWordWrap(True)
        v.addWidget(self.note)

    # ---------- Public ----------
    def set_matrix(self, A: sparse.spmatrix, spd: Optional[bool]):
        self.A = A.tocsr() if A is not None else None
        self.spd = spd
        self._refresh()

    # ---------- Helpers ----------
    def _add_section(self, title: str):
        r = self.table.rowCount()
        self.table.insertRow(r)
        item = QtWidgets.QTableWidgetItem(title)
        font = item.font(); font.setBold(True)
        item.setFont(font)
        item.setForeground(QtGui.QBrush(QtGui.QColor("#e5e7eb")))
        item.setFlags(QtCore.Qt.ItemIsEnabled)
        self.table.setItem(r, 0, item)
        # Span across all columns
        self.table.setSpan(r, 0, 1, self.table.columnCount())

    def _check_item(self, ok: Optional[bool]) -> QtWidgets.QTableWidgetItem:
        """
        ok=True  -> green ✔
        ok=False -> red ✖
        ok=None  -> gray ?
        """
        if ok is True:
            text, color = "✔", "#86efac"
        elif ok is False:
            text, color = "✖", "#fca5a5"
        else:
            text, color = "?", "#cbd5e1"
        itm = QtWidgets.QTableWidgetItem(text)
        itm.setTextAlignment(QtCore.Qt.AlignCenter)
        itm.setForeground(QtGui.QBrush(QtGui.QColor(color)))
        return itm

    def _add_row(self, name: str, typ: str, req: str, ok: Optional[bool]):
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QtWidgets.QTableWidgetItem(name))
        self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(typ))
        self.table.setItem(r, 2, QtWidgets.QTableWidgetItem(req))
        self.table.setItem(r, 3, self._check_item(ok))

    def _refresh(self):
        self.table.setRowCount(0)
        if self.A is None:
            return

        m, n = self.A.shape
        square = (m == n)
        rectangular = not square
        # cheap structural symmetry hint
        try:
            sym_hint = square and (self.A.nnz == self.A.T.nnz)
        except Exception:
            sym_hint = square
        spd = self.spd if self.spd is not None else None  # may be None

        # ---------------- Sections ----------------
        # DIRECT factorizations
        self._add_section("Direct Factorizations")
        self._add_row("Cholesky", "Direct",
                      "square, symmetric, SPD", square and (spd is True))
        self._add_row("LDLᵀ (Bunch–Kaufman)", "Direct",
                      "square, symmetric (indefinite allowed)", square and sym_hint)
        self._add_row("LU", "Direct",
                      "square (general)", square)
        self._add_row("QR (sparse/dense)", "Direct",
                      "any shape; stable; LS via QR", True)
        self._add_row("SVD (dense/iterative)", "Direct",
                      "any shape; most stable LS", True)
        self._add_row("Pivoted Cholesky (approx.)", "Direct (approx.)",
                      "SPD (approximation for kernels/low-rank)", (spd is True))

        # ITERATIVE solvers (Krylov) for Ax = b
        self._add_section("Iterative Solvers (Krylov) — Ax = b")
        self._add_row("CG (Conjugate Gradient)", "Iterative",
                      "square, symmetric, SPD", square and (spd is True))
        self._add_row("PCG (Preconditioned CG)", "Iterative",
                      "square, symmetric, SPD + preconditioner", square and (spd is True))
        self._add_row("MINRES", "Iterative",
                      "square, symmetric (SPD or indefinite)", square and sym_hint)
        self._add_row("SYMMLQ", "Iterative",
                      "square, symmetric (SPD or indefinite)", square and sym_hint)
        self._add_row("GMRES", "Iterative",
                      "square (general nonsymmetric)", square)
        self._add_row("FGMRES (flexible)", "Iterative",
                      "square (general) with variable preconditioner", square)
        self._add_row("BiCG", "Iterative",
                      "square (general); requires Aᵀ matvec", square)
        self._add_row("BiCGSTAB", "Iterative",
                      "square (general)", square)
        self._add_row("CGS", "Iterative",
                      "square (general)", square)
        self._add_row("QMR", "Iterative",
                      "square (general); often uses Aᵀ", square)
        self._add_row("TFQMR", "Iterative",
                      "square (general)", square)

        # LEAST-SQUARES solvers  (min ||Ax - b||)
        self._add_section("Least-Squares (min ||Ax − b||)")
        self._add_row("QR-based LS", "Direct",
                      "any shape; robust", True)
        self._add_row("SVD-based LS", "Direct",
                      "any shape; most robust", True)
        self._add_row("LSQR", "Iterative",
                      "any shape; works for sparse tall/skinny or wide matrices", True)
        self._add_row("LSMR", "Iterative",
                      "any shape; like LSQR, often better conditioning", True)
        self._add_row("CGNR / CGNE (aka CGLS)", "Iterative",
                      "any shape via normal equations AᵀA x = Aᵀ b (caution: conditioning)", True)

        # PRECONDITIONERS (applicability hints)
        self._add_section("Preconditioners (hints)")
        self._add_row("Jacobi / Block-Jacobi", "Preconditioner",
                      "any; cheap diagonal/block inverse", True)
        self._add_row("ILU(k)", "Preconditioner",
                      "square, general (nonsymmetric) — for GMRES/BiCGSTAB", square)
        self._add_row("IC(0) / Incomplete Cholesky", "Preconditioner",
                      "square, SPD — for CG/PCG/MINRES", square and (spd is True))
        self._add_row("Algebraic Multigrid (AMG)", "Preconditioner/Solver",
                      "best for SPD/M-matrices; symmetric problems", square and (spd is not False))
        self._add_row("Approximate Inverse (SPAI)", "Preconditioner",
                      "square (general)", square)

        # NOTES / edge cases section
        self._add_section("Notes")
        self._add_row("Normal Equations (AᵀA)", "Path",
                      "turns LS into SPD system; may worsen conditioning", True)
        self._add_row("Augmented Systems", "Path",
                      "KKT-style block solves for constrained LS/optimization", True)

        # Resize columns a bit
        self.table.resizeColumnsToContents()
