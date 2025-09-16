/******************************************************************************
 * @file constraint_set.cpp
 * @brief Assemble the sparse constraint system `C u = d` from high-level equations.
 *
 * -----------------------------------------------------------------------------
 * ## What this file does
 *
 * Takes a list of symbolic equations (each equation is a small list of
 * `(node_id, dof, coeff)` entries with a right-hand-side `rhs`) and a global
 * DOF map, and turns them into a concrete sparse matrix/vector pair:
 *
 *      C ∈ R^{m×n},   d ∈ R^m   such that   C u = d .
 *
 * Each equation produces one row in `C` (unless dropped by the zero-row policy).
 *
 * -----------------------------------------------------------------------------
 * ## ELI5
 *
 * Think of each equation as “put these numbers into these columns, all on the
 * same row, and set the row’s right-hand side to this value”. The DOF map tells
 * you *which column* corresponds to `(node_id, dof)`. If a DOF is inactive
 * (`-1`), we simply skip it.
 *
 * -----------------------------------------------------------------------------
 * ## Scaling (robustness)
 *
 * - **Column scaling** (default on): scale each nonzero column of `C` to unit
 *   2-norm. This improves the numerical behavior of the QR factorization used
 *   later to build the null-space map (ConstraintBuilder).
 *
 * - **Row scaling** (default off): optional left scaling to unit ∞-norm to
 *   compress the dynamic range of each row.
 *
 * The applied scale factors are stored in `col_scale` and `row_scale` so you
 * can track/report or (if ever needed) undo the scaling externally.
 *
 * -----------------------------------------------------------------------------
 * ## Row dropping
 *
 * If an equation touches only inactive DOFs, it would assemble to a **zero row**.
 * By default (zero_row_drop_tol == 0), we **keep** it (useful for diagnostics).
 * If `zero_row_drop_tol > 0`, such zero rows are **dropped** entirely.
 *
 * -----------------------------------------------------------------------------
 * @see constraint_set.h  (API and high-level documentation)
 * @see equation.h        (symbolic equation representation)
 * @date    14.09.2025
 * @author  Finn
 ******************************************************************************/

#include "constraint_set.h"
#include "../core/logging.h"
#include <unordered_map>

namespace fem::constraint {

void ConstraintSet::assemble(const SystemDofIds& dofs, Index n_dofs)
{
    // Keep a reference to the DOF map and declare number of columns (total DOFs).
    dof_map = &dofs;
    n = n_dofs;

    // Triplet buffer for sparse C; heuristic reserve keeps reallocation low
    TripletList trips;
    trips.reserve(64 * equations.size());

    // Bookkeeping (which original equation produced row i in C)
    kept_row_ids.clear();
    kept_row_ids.reserve(equations.size());

    // Temporary RHS buffer (we’ll map it into Eigen later)
    std::vector<Precision> d_vals;
    d_vals.reserve(equations.size());

    // --------------------- Pass 1: build rows --------------------------------
    Index row = 0;
    for (Index ei = 0; ei < static_cast<Index>(equations.size()); ++ei) {
        const auto& eq = equations[ei];
        bool has_entry = false;

        // For each entry in the equation, find its global column and add a triplet
        for (const auto& t : eq.entries) {
            const int dof_id = dofs(t.node_id, t.dof);
            if (dof_id < 0) continue;                // inactive DOF → skip
            has_entry = true;
            trips.emplace_back(row, dof_id, t.coeff);
        }

        // Decide whether to keep the row:
        //  - If it has at least one active entry, we keep it.
        //  - If it is empty and zero_row_drop_tol == 0, we still keep it as a zero row.
        //  - If it is empty and zero_row_drop_tol > 0, we drop it.
        if (has_entry || opt.zero_row_drop_tol <= Precision(0)) {
            kept_row_ids.push_back(ei);
            d_vals.push_back(eq.rhs);
            ++row;
        }
        // else: dropped, do not increment row (no matrix entries added for it)
    }

    m = row;

    // Edge case: nothing assembled
    if (m == 0) {
        C.resize(0, n);
        d.resize(0);
        col_scale = DynamicVector::Ones(n);
        row_scale = DynamicVector::Ones(m);
        logging::info(true, "[ConstraintSet] No effective constraint rows.");
        return;
    }

    // --------------------- Build C and d -------------------------------------
    C.resize(m, n);
    C.setFromTriplets(trips.begin(), trips.end());
    d = Eigen::Map<DynamicVector>(d_vals.data(), static_cast<Eigen::Index>(d_vals.size()));

    // Initialize scaling vectors to identity
    col_scale = DynamicVector::Ones(n);
    row_scale = DynamicVector::Ones(m);

    // --------------------- Column scaling (2-norm) ----------------------------
    if (opt.scale_columns) {
        for (int j = 0; j < C.outerSize(); ++j) {
            Precision norm2 = 0;
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                norm2 += it.value() * it.value();
            }
            const Precision s = (norm2 > 0) ? Precision(1.0) / std::sqrt(norm2) : Precision(1.0);
            if (s != Precision(1.0)) {
                for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                    it.valueRef() *= s;
                }
            }
            col_scale[j] = s;
        }
    }

    // --------------------- Row scaling (∞-norm) -------------------------------
    if (opt.scale_rows) {
        std::vector<Precision> row_inf(m, 0);
        // Compute max absolute entry per row
        for (int j = 0; j < C.outerSize(); ++j) {
            for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
                row_inf[it.row()] = std::max(row_inf[it.row()], std::abs(it.value()));
            }
        }

        // Build diagonal scaling D and apply: C ← D C, d ← D d
        SparseMatrix D(m, m);
        D.reserve(Eigen::VectorXi::Constant(m, 1));
        for (int i = 0; i < m; ++i) {
            const Precision s = (row_inf[i] > 0) ? Precision(1.0) / row_inf[i] : Precision(1.0);
            D.insert(i, i) = s;
            row_scale[i]   = s;
        }
        D.makeCompressed();

        C = D * C;
        d = row_scale.asDiagonal() * d;
    }

//    // --------------------- Logging -------------------------------------------
//    logging::info(true,
//        "[ConstraintSet] Assembled C (", m, "x", n, ", nnz=", C.nonZeros(),
//        "), d (m=", m, "), kept_rows=", kept_row_ids.size());
}

} // namespace fem::constraint
