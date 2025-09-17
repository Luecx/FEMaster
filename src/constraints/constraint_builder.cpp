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
#include "constraint_map.h"
#include "../core/logging.h"

#include <Eigen/SparseQR>                        // SparseQR
#include <Eigen/OrderingMethods>                 // COLAMDOrdering
#include <algorithm>
#include <unordered_map>
#include <cmath>

namespace fem::constraint {

/*** Helpers ******************************************************************/

static void compress_zero_columns(const SparseMatrix& C,
                                  std::vector<int>&   used_cols,
                                  Eigen::VectorXi&    old2new,
                                  SparseMatrix&       C_use)
{
    const int n = C.cols();
    used_cols.clear();
    used_cols.reserve(n);

    for (int j = 0; j < C.outerSize(); ++j) {
        bool nonzero = false;
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) { nonzero = true; break; }
        if (nonzero) used_cols.push_back(j);
    }

    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < (int)used_cols.size(); ++k) old2new[used_cols[k]] = k;

    C_use.resize(C.rows(), (int)used_cols.size());
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C.nonZeros());
    for (int jnew = 0; jnew < (int)used_cols.size(); ++jnew) {
        const int j = used_cols[jnew];
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            trips.emplace_back(it.row(), jnew, it.value());
        }
    }
    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

/*** Builder ******************************************************************/

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt;
    return build(set, opt);
}

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set, const Options& opt)
{
    using Mat  = Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec  = Eigen::Matrix<Precision, Eigen::Dynamic, 1>;
    using Perm = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;

    // ---- Compress to constrained DOFs only ---------------------------------
    std::vector<int> used;
    Eigen::VectorXi old2new;
    SparseMatrix C_use;
    compress_zero_columns(set.C, used, old2new, C_use);

    ConstraintMap M;
    Report rep;
    rep.m = set.m;
    rep.n = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    const int n_use = (int)used.size();

    // Trivial (no constrained columns): identity map
    if (n_use == 0) {
        M.n_  = set.n;
        M.r_  = 0;
        M.nm_ = set.n;
        M.masters_.resize(set.n);
        M.slaves_.clear();
        for (int j = 0; j < set.n; ++j) M.masters_[j] = j;

        M.T_ = SparseMatrix(set.n, set.n); M.T_.setIdentity();
        M.X_.resize(0, set.n);
        M.u_p_ = DynamicVector::Zero(set.n);

        rep.rank = 0;
        rep.n_redundant_rows = rep.m;
        rep.feasible = true;
        rep.d_norm = (set.d.size() == 0) ? 0 : set.d.norm();
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
    Perm P = sqr.colsPermutation();
    SparseMatrix Rsp = sqr.matrixR();   // sparse upper-triangular (m x n_use)

    // Rank detection via |diag(R)| relative to max |diag|
    const int rmax = std::min<int>(Rsp.rows(), Rsp.cols());
    Precision max_abs_diag = 0;
    for (int k = 0; k < rmax; ++k) {
        const Precision ad = std::abs(Rsp.coeff(k, k)); // diagonal exists for QR
        if (ad > max_abs_diag) max_abs_diag = ad;
    }
    rep.R11_max_diag = max_abs_diag;

    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);
    int r = 0;
    for (; r < rmax; ++r) {
        if (std::abs(Rsp.coeff(r, r)) <= tol) break;
    }
    rep.rank = r;
    const int nm_use = n_use - r;

    // ---- Build Xd = -R11^{-1} R12  (thin dense, r x nm_use) -----------------
    Mat Xd = Mat::Zero(r, nm_use);
    if (r > 0 && nm_use > 0) {
        // Extract R11 and R12 as small dense from sparse R
        Mat R11 = Mat::Zero(r, r);
        Mat R12 = Mat::Zero(r, nm_use);

        // Iterate by column j of Rsp (upper-triangular)
        const int Rcols = Rsp.cols();
        for (int j = 0; j < Rcols; ++j) {
            for (SparseMatrix::InnerIterator it(Rsp, j); it; ++it) {
                const int i = it.row();           // row index
                if (i >= r) continue;             // only rows [0..r-1] matter
                const Precision v = it.value();
                if (j < r) {
                    R11(i, j) = v;
                } else if (j < r + nm_use) {
                    R12(i, j - r) = v;
                } // else: beyond the leading [r | nm_use] block; ignore
            }
        }

        // Solve X = -R11^{-1} R12 (upper-triangular solve; no manual regularization)
        Xd = (-R11.template triangularView<Eigen::Upper>().solve(R12)).eval();
    }

    // ---- Partition columns by permutation P --------------------------------
    const auto& p = P.indices(); // P*e_i = e_{p(i)}

    std::vector<int> slaves_loc;  slaves_loc.reserve(r);
    std::vector<int> masters_loc; masters_loc.reserve(nm_use);
    for (int i = 0; i < r;      ++i) slaves_loc.push_back(  p(i)     );
    for (int j = 0; j < nm_use; ++j) masters_loc.push_back( p(r + j) );

    // ---- Build T and X directly (no dense Zp/Tuse) --------------------------
    std::vector<char> is_used(set.n, 0);
    for (int c : used) is_used[c] = 1;
    int n_never = 0; for (int g = 0; g < set.n; ++g) if (!is_used[g]) ++n_never;
    const int nm_total = nm_use + n_never;

    // global index lists
    std::vector<Index> masters_glob; masters_glob.reserve(nm_total);
    std::vector<Index> slaves_glob;  slaves_glob.reserve(r);
    for (int j = 0; j < nm_use; ++j) masters_glob.push_back( used[ masters_loc[j] ] );
    for (int g = 0; g < set.n; ++g) if (!is_used[g]) masters_glob.push_back(g);
    for (int i = 0; i < r;      ++i) slaves_glob.push_back( used[ slaves_loc[i] ] );

    // T_full (n x nm_total)
    SparseMatrix T_full(set.n, nm_total);
    {
        std::vector<Eigen::Triplet<Precision>> tTrips;
        // Rough pre-reserve: one 1 per constrained master + average sparsity for X + identities for never-used
        tTrips.reserve( (size_t)nm_use * (size_t)(1 + std::max(1, r/16)) + (size_t)n_never );

        // constrained master columns
        for (int j = 0; j < nm_use; ++j) {
            const int row_master_loc = masters_loc[j];
            tTrips.emplace_back( used[row_master_loc], j, Precision(1) );
            for (int i = 0; i < r; ++i) {
                Precision v = Xd(i, j);
                if (v != Precision(0)) tTrips.emplace_back( used[ slaves_loc[i] ], j, v );
            }
        }
        // append identity for never-used masters
        int col_off = nm_use;
        for (int g = 0; g < set.n; ++g)
            if (!is_used[g]) tTrips.emplace_back(g, col_off++, Precision(1));

        T_full.setFromTriplets(tTrips.begin(), tTrips.end());
        T_full.makeCompressed();
    }

    // X_ (r x nm_total)
    SparseMatrix X(r, nm_total);
    if (r > 0 && nm_use > 0) {
        std::vector<Eigen::Triplet<Precision>> xTrips;
        xTrips.reserve( (size_t)r * (size_t)nm_use / 4 + 1 ); // heuristic
        for (int j = 0; j < nm_use; ++j)
            for (int i = 0; i < r; ++i) {
                Precision v = Xd(i, j);
                if (v != Precision(0)) xTrips.emplace_back(i, j, v);
            }
        X.setFromTriplets(xTrips.begin(), xTrips.end());
        X.makeCompressed();
    } else {
        X.resize(r, nm_total); // empty
    }

    // ---- Populate map -------------------------------------------------------
    M.n_       = set.n;
    M.r_       = r;
    M.nm_      = (Index)masters_glob.size();
    M.masters_ = std::move(masters_glob);
    M.slaves_  = std::move(slaves_glob);
    M.T_       = std::move(T_full);
    M.X_       = std::move(X);
    M.u_p_     = DynamicVector::Zero(set.n);

    rep.feasible          = true;
    rep.d_norm            = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows  = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

#ifndef NDEBUG
    // ---- Invariants (debug only) --------------------------------------------
    // 1) C*T ≈ 0
    {
        const int nm = (int)M.nm_;
        DynamicVector e = DynamicVector::Zero(nm);
        double max_norm = 0.0;
        for (int k = 0; k < nm; ++k) {
            e.setZero(); e[k] = 1;
            DynamicVector cu = set.C * (M.T_ * e);
            max_norm = std::max(max_norm, (double)cu.norm());
        }
        logging::error(max_norm <= 1e-12, "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ", max_norm);
    }
    // 2) apply_T vs explicit T_
    {
        DynamicVector q = DynamicVector::Random(M.nm_);
        DynamicVector u_explicit = M.T_ * q + M.u_p_;
        DynamicVector u_fast     = DynamicVector::Zero(M.n_);
        M.apply_T(q, u_fast);
        logging::error( (u_explicit - u_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_T mismatch explicit T.");
    }
    // 3) apply_Tt vs explicit T_^T
    {
        DynamicVector y = DynamicVector::Random(M.n_);
        DynamicVector z_explicit = M.T_.transpose() * y;
        DynamicVector z_fast     = DynamicVector::Zero(M.nm_);
        M.apply_Tt(y, z_fast);
        logging::error( (z_explicit - z_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    }
#endif

    // ---- Inhomogeneous: particular solution u_p (min-norm) ------------------
    if (!rep.homogeneous) {
        Vec d_dense(set.m); for (int i = 0; i < set.m; ++i) d_dense[i] = set.d[i];

        // Solve min-norm u_use for C_use * u_use ≈ d
        Vec u_use = sqr.solve(d_dense);

        // Lift to full vector (zeros elsewhere)
        DynamicVector u_any = DynamicVector::Zero(set.n);
        for (int k = 0; k < n_use; ++k) u_any[used[k]] = u_use[k];

        // Residual in original scaling
        DynamicVector resid = set.C * u_any - set.d;
        rep.residual_norm   = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            // Project u_any onto affine map to get u_p (store only slave offsets)
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any;
            for (int i = 0; i < M.r_; ++i) {
                const Index gi = M.slaves_[i];
                M.u_p_[gi] = u_any[gi] - Xq[i];
            }
        }
    }

    return {M, rep};
}

} // namespace fem::constraint
