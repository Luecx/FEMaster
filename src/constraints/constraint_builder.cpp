/******************************************************************************
 * @file constraint_builder.cpp
 * @brief Build a reduced null-space map u = u_p + T q from C u = d.
 *
 * Sparse pipeline (no dense QR, no dense intermediates):
 *   - Compress zero columns of C → C_use (m x n_use)
 *   - SparseQR (COLAMD): C_use * P = Q * R
 *   - Rank r from diag(R) via relative threshold
 *   - Solve X = -R11^{-1} R12    (X held as thin dense r×(n_use−r))
 *   - Build T and X directly from (used, P, X) as sparse triplets
 *   - Inhomogeneous u_p via SparseQR::solve (min-norm), then slave offsets
 *
 * All application/assembly paths remain sparse-only (see ConstraintMap).
 *
 * @date    14.09.2025  (SparseQR-only fast path)
 * @author  Finn
 ******************************************************************************/

#include "constraint_builder.h"

#include "../core/logging.h"
#include "constraint_map.h"

#include <Eigen/OrderingMethods>    // COLAMDOrdering
#include <Eigen/SparseQR>           // SparseQR
#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace fem::constraint {

/*** Helpers ******************************************************************/

static void compress_zero_columns(const SparseMatrix& C,
                                  std::vector<int>&   used_cols,
                                  Eigen::VectorXi&    old2new,
                                  SparseMatrix&       C_use) {
    const int n = C.cols();
    used_cols.clear();
    used_cols.reserve(n);

    for (int j = 0; j < C.outerSize(); ++j) {
        bool nonzero = false;
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            nonzero = true;
            break;
        }
        if (nonzero)
            used_cols.push_back(j);
    }

    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < (int) used_cols.size(); ++k)
        old2new[used_cols[k]] = k;

    C_use.resize(C.rows(), (int) used_cols.size());
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C.nonZeros());
    for (int jnew = 0; jnew < (int) used_cols.size(); ++jnew) {
        const int j = used_cols[jnew];
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            trips.emplace_back(it.row(), jnew, it.value());
        }
    }
    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

/*** Builder ******************************************************************/

std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt;
    return build(set, opt);
}

std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set,
                                                                             const Options&       opt) {
    using Mat  = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec  = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;
    using Perm = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;

    // ---- Compress to constrained DOFs only ---------------------------------
    std::vector<int> used;
    Eigen::VectorXi  old2new;
    SparseMatrix     C_use;
    compress_zero_columns(set.C, used, old2new, C_use);

    ConstraintMap M;
    Report        rep;
    rep.m           = set.m;
    rep.n           = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    const int n_use = (int) used.size();

    // Trivial (no constrained columns): identity map
    if (n_use == 0) {
        M.n_  = set.n;
        M.r_  = 0;
        M.nm_ = set.n;
        M.masters_.resize(set.n);
        M.slaves_.clear();
        for (int j = 0; j < set.n; ++j)
            M.masters_[j] = j;

        M.T_ = SparseMatrix(set.n, set.n);
        M.T_.setIdentity();
        M.X_.resize(0, set.n);
        M.u_p_               = DynamicVector::Zero(set.n);

        rep.rank             = 0;
        rep.n_redundant_rows = rep.m;
        rep.feasible         = true;
        rep.d_norm           = (set.d.size() == 0) ? 0 : set.d.norm();
        return {M, rep};
    }

    // ---- Sparse, rank-revealing QR:  C_use * P = Q * R ---------------------
    // (optional) clean up zeros to help QR a bit
    C_use.prune(0.0);
    C_use.makeCompressed();
    // AMD is typically faster here than COLAMD; pivoting fully off.
    // (We determine rank ourselves via |diag(R)|.)
    Eigen::SparseQR<SparseMatrix, Eigen::AMDOrdering<int>> sqr;
    sqr.setPivotThreshold(0.0);
    sqr.compute(C_use);

    logging::error(sqr.info() == Eigen::Success, "[ConstraintBuilder] SparseQR failed.");

    // Extract permutation and R
    Perm         P   = sqr.colsPermutation();
    SparseMatrix Rsp = sqr.matrixR();    // sparse upper-triangular (m x n_use)

    // Rank detection via |diag(R)| relative to max |diag|
    const int rmax         = std::min<int>(Rsp.rows(), Rsp.cols());
    Precision max_abs_diag = 0;
    for (int k = 0; k < rmax; ++k) {
        const Precision ad = std::abs(Rsp.coeff(k, k));    // diagonal exists for QR
        if (ad > max_abs_diag)
            max_abs_diag = ad;
    }
    rep.R11_max_diag    = max_abs_diag;

    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);
    int             r   = 0;
    for (; r < rmax; ++r) {
        if (std::abs(Rsp.coeff(r, r)) <= tol)
            break;
    }
    rep.rank         = r;
    const int nm_use = n_use - r;

    // ---- Partition columns by permutation P --------------------------------
    const auto&      p = P.indices();    // P*e_i = e_{p(i)}
    std::vector<int> slaves_loc;
    slaves_loc.reserve(r);
    std::vector<int> masters_loc;
    masters_loc.reserve(nm_use);
    for (int i = 0; i < r; ++i)
        slaves_loc.push_back(p(i));
    for (int j = 0; j < nm_use; ++j)
        masters_loc.push_back(p(r + j));

    // ---- Build T and X using sparse triangular solves (no dense Xd) --------
    std::vector<char> is_used(set.n, 0);
    for (int c : used)
        is_used[c] = 1;
    int n_never = 0;
    for (int g = 0; g < set.n; ++g)
        if (!is_used[g])
            ++n_never;
    const int nm_total = nm_use + n_never;

    // global index lists
    std::vector<Index> masters_glob;
    masters_glob.reserve(nm_total);
    std::vector<Index> slaves_glob;
    slaves_glob.reserve(r);
    for (int j = 0; j < nm_use; ++j)
        masters_glob.push_back(used[masters_loc[j]]);
    for (int g = 0; g < set.n; ++g)
        if (!is_used[g])
            masters_glob.push_back(g);
    for (int i = 0; i < r; ++i)
        slaves_glob.push_back(used[slaves_loc[i]]);

    Rsp.makeCompressed();

    // Build sparse R11 (upper) and access R12 via columns r..r+nm_use-1
    SparseMatrix R11s(r, r);
    {
        std::vector<Eigen::Triplet<Precision>> trips;
        trips.reserve((size_t) r * 4);    // heuristic
        for (int j = 0; j < r; ++j) {
            for (SparseMatrix::InnerIterator it(Rsp, j); it; ++it) {
                const int i = it.row();
                if (i > j)
                    break;    // upper
                trips.emplace_back(i, j, it.value());
            }
        }
        R11s.setFromTriplets(trips.begin(), trips.end());
        R11s.makeCompressed();
    }

    // Triplets for T and X we will fill streaming
    std::vector<Eigen::Triplet<Precision>> tTrips;
    std::vector<Eigen::Triplet<Precision>> xTrips;
    // Rough reserve: identity on nm_use + few slave coeffs + identities for never-used
    tTrips.reserve((size_t) nm_use + (size_t) std::max(1, r / 16) * (size_t) nm_use + (size_t) n_never);
    xTrips.reserve((size_t) std::max(1, r / 16) * (size_t) nm_use);

    // Work vectors for solve
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> b(r), x(r);

    // Per master column j: solve R11 * x = -R12(:,j) sparsely and emit trips
    for (int j = 0; j < nm_use; ++j) {
        b.setZero();

        // RHS b = -R12(:,j) lives in column (r + j) of R
        const int colR = r + j;
        if (colR < Rsp.cols()) {
            for (SparseMatrix::InnerIterator it(Rsp, colR); it; ++it) {
                const int i = it.row();
                if (i >= r)
                    break;    // only rows in top block
                b[i] = -it.value();
            }
        }

        // Sparse triangular solve: R11 (upper) * x = b
        x = R11s.template triangularView<Eigen::Upper>().solve(b);

        // Optional per-column drop (relative to column norm)
        if (opt.threshold_X) {
            const Precision cn = x.norm();
            if (cn > Precision(0)) {
                const Precision thr = cn * opt.X_drop_tol_rel;
                for (int i = 0; i < r; ++i)
                    if (std::abs(x[i]) < thr)
                        x[i] = Precision(0);
            }
        }

        // Global master column index (first nm_use columns are constrained masters)
        const int jglob = j;

        // Identity on the master row of this column
        const int row_master_loc = masters_loc[j];
        tTrips.emplace_back(used[row_master_loc], jglob, Precision(1));

        // Add slave couplings for this column
        for (int i = 0; i < r; ++i) {
            const Precision v = x[i];
            if (v == Precision(0))
                continue;
            const int row_slave_use = slaves_loc[i];
            tTrips.emplace_back(used[row_slave_use], jglob, v);
            // X row i aligns with slaves_glob[i], column jglob matches T's
            xTrips.emplace_back(i, jglob, v);
        }
    }

    // Append identity for never-used masters
    {
        int col_off = nm_use;
        for (int g = 0; g < set.n; ++g)
            if (!is_used[g])
                tTrips.emplace_back(g, col_off++, Precision(1));
    }

    // Finalize T and X
    SparseMatrix T_full(set.n, nm_total);
    T_full.setFromTriplets(tTrips.begin(), tTrips.end());
    T_full.makeCompressed();

    SparseMatrix X(r, nm_total);
    X.setFromTriplets(xTrips.begin(), xTrips.end());
    X.makeCompressed();

    // ---- Populate map -------------------------------------------------------
    M.n_                 = set.n;
    M.r_                 = r;
    M.nm_                = (Index) masters_glob.size();
    M.masters_           = std::move(masters_glob);
    M.slaves_            = std::move(slaves_glob);
    M.T_                 = std::move(T_full);
    M.X_                 = std::move(X);
    M.u_p_               = DynamicVector::Zero(set.n);

    rep.feasible         = true;
    rep.d_norm           = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

#ifndef NDEBUG
    // ---- Invariants (debug only) --------------------------------------------
    // 1) C*T ≈ 0
    {
        const int     nm       = (int) M.nm_;
        DynamicVector e        = DynamicVector::Zero(nm);
        double        max_norm = 0.0;
        for (int k = 0; k < nm; ++k) {
            e.setZero();
            e[k]             = 1;
            DynamicVector cu = set.C * (M.T_ * e);
            max_norm         = std::max(max_norm, (double) cu.norm());
        }
        logging::error(max_norm <= 1e-12, "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ", max_norm);
    }
    // 2) apply_T vs explicit T_
    {
        DynamicVector q          = DynamicVector::Random(M.nm_);
        DynamicVector u_explicit = M.T_ * q + M.u_p_;
        DynamicVector u_fast     = DynamicVector::Zero(M.n_);
        M.apply_T(q, u_fast);
        logging::error((u_explicit - u_fast).norm() <= 1e-12, "[ConstraintBuilder] apply_T mismatch explicit T.");
    }
    // 3) apply_Tt vs explicit T_^T
    {
        DynamicVector y          = DynamicVector::Random(M.n_);
        DynamicVector z_explicit = M.T_.transpose() * y;
        DynamicVector z_fast     = DynamicVector::Zero(M.nm_);
        M.apply_Tt(y, z_fast);
        logging::error((z_explicit - z_fast).norm() <= 1e-12, "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    }
#endif

    // ---- Inhomogeneous: particular solution u_p (min-norm) ------------------
    if (!rep.homogeneous) {
        Vec d_dense(set.m);
        for (int i = 0; i < set.m; ++i)
            d_dense[i] = set.d[i];

        // Solve min-norm u_use for C_use * u_use ≈ d
        Vec u_use = sqr.solve(d_dense);

        // Lift to full vector (zeros elsewhere)
        DynamicVector u_any = DynamicVector::Zero(set.n);
        for (int k = 0; k < n_use; ++k)
            u_any[used[k]] = u_use[k];

        // Residual in original scaling
        DynamicVector resid      = set.C * u_any - set.d;
        rep.residual_norm        = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible             = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            // Project u_any onto affine map to get u_p (store only slave offsets)
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j)
                q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any;
            for (int i = 0; i < M.r_; ++i) {
                const Index gi = M.slaves_[i];
                M.u_p_[gi]     = u_any[gi] - Xq[i];
            }
        }
    }

    return {M, rep};
}

}    // namespace fem::constraint
