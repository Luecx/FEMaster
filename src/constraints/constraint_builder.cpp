/******************************************************************************
 * @file constraint_builder.cpp
 * @brief Implementation: build a reduced null-space map u = u_p + T q from C u = d.
 *
 * This file contains the robust, rank-revealing construction used to transform
 * linear constraints
 *
 *     C u = d,   with  C ∈ R^{m×n}, d ∈ R^m
 *
 * into an explicit **affine map** of the form
 *
 *     u = u_p + T q,
 *
 * where:
 *   • q ∈ R^{n−r} are the **independent (master) DOFs**,
 *   • the remaining r DOFs are **slaves** expressed via X (internal to T),
 *   • u_p is a **particular solution** if d ≠ 0 (offset in the constraint manifold).
 *
 * Once built, the map allows you to reduce systems like K u = f to
 *
 *     A q = b,  with  A = Tᵀ K T,  b = Tᵀ (f − K u_p),
 *
 * solve for q, and recover u = u_p + T q that **exactly** satisfies C u = d
 * (for feasible d, up to floating point tolerance).
 *
 * -----------------------------------------------------------------------------
 * ## High-level algorithm (ELI5 meets math)
 *
 * 1) Factorize C using rank-revealing QR with column pivoting:
 *
 *        C P = Q R
 *
 *    P reorders columns (DOFs) so that the well-determined DOFs appear first.
 *    R is (upper) triangular; its first r diagonals above a tolerance define
 *    the **numerical rank r** = rank(C).
 *
 * 2) Split the columns into **slaves** and **masters** by that permutation:
 *        [u_S; u_M] := Pᵀ u   with  |u_S| = r, |u_M| = n − r.
 *
 * 3) Read the top r rows of R:  R11 u_S + R12 u_M = d′  (with d′ = Qᵀ d).
 *    For homogeneous constraints (d = 0):  u_S = −R11⁻¹ R12 u_M.
 *    For inhomogeneous constraints (d ≠ 0): also build a particular solution u_p.
 *
 * 4) Convert back to original DOF ordering and assemble:
 *        X  := −R11⁻¹ R12,
 *        T  := [ I (on masters’ rows) ; X (scattered on slaves’ rows) ],
 *        u_p (zeros for homogeneous).
 *
 * 5) Provide a **Report** (rank, feasibility, suspect rows) alongside a
 *    **ConstraintMap** that owns T, X, u_p and efficient application helpers.
 *
 * -----------------------------------------------------------------------------
 * ## Robustness notes
 *
 * • We use a **dense fallback** (ColPivHouseholderQR) when m is small, which
 *   tends to be more reliable than sparse QR for tiny systems and is cheap.
 * • Rank is decided by comparing |R_kk| to `rank_tol_rel * max_k |R_kk|`.
 * • Infeasibility (contradictory constraints) is detected for d ≠ 0 by
 *   measuring the least-squares residual ‖C u_any − d‖ and listing suspect rows.
 *
 * @see constraint_builder.h for a friendly overview and usage example.
 ******************************************************************************/

#include "constraint_builder.h"
#include "constraint_map.h"
#include "../core/logging.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <Eigen/QR>   // ColPivHouseholderQR
#include <algorithm>
#include <sstream>

namespace fem::constraint {

using QR = Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>;

/******************************************************************************
 * @brief Read |R(k,k)| from a sparse, column-major, upper-triangular block.
 *
 * This utility iterates the k-th column of R and returns the absolute value of
 * the diagonal entry. If absent (should not happen for valid triangular data),
 * it returns 0.
 *
 * @param R Column-major sparse matrix that stores an upper-triangular factor.
 * @param k Diagonal index.
 * @return |R(k,k)| as a Precision.
 ******************************************************************************/
static Precision read_diag_abs(const SparseMatrix& R, int k) {
    for (SparseMatrix::InnerIterator it(R, k); it; ++it) {
        if (it.row() == k) return std::abs(it.value());
        if (it.row() >  k) break; // passed the diagonal in this column
    }
    return Precision(0);
}

/******************************************************************************
 * @brief Assemble (X, T) and index partitions from R and the column permutation.
 *
 * Assumes the top `r` rows of R have already been extracted and passed in.
 * The column permutation `permCols` maps **permuted** column index -> original DOF.
 *
 * Partition:
 *   • first r permuted columns → **slaves** (size r),
 *   • remaining n − r columns → **masters** (size nm).
 *
 * Then solve X = −R11⁻¹ R12 by upper-triangular back-substitution column-wise,
 * and scatter (masters → identity rows, slaves ← X rows) to build T.
 *
 * @param n        Total number of DOFs (columns of C).
 * @param r        Numerical rank.
 * @param permCols Column permutation from QR (permuted col -> original col).
 * @param R        Top r rows of the R-factor (r × n) as a sparse matrix.
 * @param slaves   (out) global DOF indices marked as slaves, size r.
 * @param masters  (out) global DOF indices marked as masters, size n − r.
 * @param X        (out) r × (n − r) sparse matrix with slave coefficients.
 * @param T        (out) n × (n − r) sparse map from masters to full vector.
 ******************************************************************************/
static void build_T_from_R_and_perm(
    int n, int r, const Eigen::VectorXi& permCols, const SparseMatrix& R,
    std::vector<Index>& slaves, std::vector<Index>& masters,
    SparseMatrix& X, SparseMatrix& T)
{
    const int nm = n - r;

    // 1) Column partition via permutation: first r → slaves, rest → masters.
    slaves.resize(r);
    masters.resize(nm);
    for (int k = 0; k < r; ++k) slaves[k] = permCols[k];
    for (int k = r; k < n; ++k) masters[k - r] = permCols[k];

    // 2) Extract R11 (r×r) and R12 (r×nm).
    SparseMatrix R11 = R.leftCols(r);
    SparseMatrix R12 = R.rightCols(nm);

    // 3) Solve X = −R11⁻¹ R12 (upper-triangular). Build X by columns.
    X.resize(r, nm);
    std::vector<Eigen::Triplet<Precision>> xTrips;
    xTrips.reserve(std::max<int>(1, R12.nonZeros()));

    for (int j = 0; j < R12.outerSize(); ++j) {
        DynamicVector b = DynamicVector::Zero(r);
        for (SparseMatrix::InnerIterator it(R12, j); it; ++it) b[it.row()] = -it.value();
        DynamicVector x = (r > 0)
            ? R11.template triangularView<Eigen::Upper>().solve(b)
            : DynamicVector();
        for (int i = 0; i < x.size(); ++i) {
            if (x[i] != Precision(0)) xTrips.emplace_back(i, j, x[i]);
        }
    }
    X.setFromTriplets(xTrips.begin(), xTrips.end());
    X.makeCompressed();

    // 4) Assemble T: identity on master rows; X on slave rows.
    T.resize(n, nm);
    std::vector<Eigen::Triplet<Precision>> tTrips;
    tTrips.reserve(std::max<int>(1, nm) + X.nonZeros());

    // Masters: put the identity rows.
    for (int j = 0; j < nm; ++j) tTrips.emplace_back(masters[j], j, Precision(1));

    // Slaves: scatter X rows to the correct original DOF rows.
    for (int j = 0; j < X.outerSize(); ++j) {
        for (SparseMatrix::InnerIterator it(X, j); it; ++it) {
            const Index i = it.row();                // row in X = index in slave list
            tTrips.emplace_back(slaves[i], j, it.value());
        }
    }

    T.setFromTriplets(tTrips.begin(), tTrips.end());
    T.makeCompressed();
}

/******************************************************************************
 * @brief Dense QR fallback for small m: robust rank and clean X/T/u_p assembly.
 *
 * For small numbers of constraints (tiny m), dense ColPivHouseholderQR provides
 * very reliable pivoting and rank decisions and is cheap. We keep an entirely
 * consistent code path by converting the top r rows of the (dense) R to a
 * sparse representation before calling the common T-builder.
 *
 * @param set Assembled constraint set (contains C, d, sizes).
 * @param opt Numerical options (rank tolerance, feasibility tolerance, etc.).
 * @return {ConstraintMap, Report}
 ******************************************************************************/
static std::pair<ConstraintMap, ConstraintBuilder::Report>
build_dense_qr(const ConstraintSet& set, const ConstraintBuilder::Options& opt)
{
    ConstraintMap M;
    ConstraintBuilder::Report rep;

    rep.m = set.m;
    rep.n = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    // 1) Factorize C densely with column pivoting.
    Eigen::MatrixXd Cd = Eigen::MatrixXd(set.C);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dqr(Cd);

    // Extract an upper-triangular view of matrixR().
    Eigen::MatrixXd R = dqr.matrixR()
                            .topLeftCorner(set.m, set.n)
                            .template triangularView<Eigen::Upper>();

    // 2) Convert top m rows of R to a sparse matrix Rs for uniform downstream use.
    SparseMatrix Rs(set.m, set.n);
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(set.C.nonZeros());
    for (int i = 0; i < R.rows(); ++i) {
        for (int j = i; j < R.cols(); ++j) {
            double v = R(i, j);
            if (v != 0.0) trips.emplace_back(i, j, v);
        }
    }
    Rs.setFromTriplets(trips.begin(), trips.end());
    Rs.makeCompressed();

    // 3) Rank decision from diag(R) with robust relative tolerance.
    Precision max_abs_diag = 0;
    std::vector<Precision> dabs(std::min(set.m, set.n), 0);
    for (int k = 0; k < static_cast<int>(dabs.size()); ++k) {
        dabs[k] = std::abs(R(k, k));
        if (dabs[k] > max_abs_diag) max_abs_diag = dabs[k];
    }
    rep.R11_max_diag = max_abs_diag;
    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);

    int r = 0; // count all diagonals above tol (do NOT early-break)
    for (int k = 0; k < static_cast<int>(dabs.size()); ++k) {
        if (dabs[k] > tol) ++r;
    }
    rep.rank = r;

    // 4) Column permutation and build of (X, T) and partition indices.
    Eigen::VectorXi p = dqr.colsPermutation().indices();
    build_T_from_R_and_perm(set.n, r, p, Rs.topRows(r), M.slaves_, M.masters_, M.X_, M.T_);

    // 5) Populate sizes & defaults; assess feasibility if d ≠ 0.
    M.n_   = set.n;
    M.r_   = r;
    M.nm_  = set.n - r;
    M.u_p_ = DynamicVector::Zero(set.n);
    rep.feasible = true;
    rep.d_norm   = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

    if (!rep.homogeneous) {
        // Least-squares particular solution (dense).
        DynamicVector u_any = dqr.solve(set.d);
        DynamicVector resid = set.C * u_any - set.d;
        rep.residual_norm   = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            // Build u_p so that any constrained solution can be written as u = u_p + T q.
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any;      // in slave space
            DynamicVector s  = DynamicVector::Zero(M.r_);
            for (int i = 0; i < M.r_; ++i) s[i] = u_any[M.slaves_[i]] - Xq[i];
            for (int i = 0; i < M.r_; ++i) M.u_p_[M.slaves_[i]] = s[i];

//            logging::info(true, "[ConstraintBuilder] Built particular solution u_p (dense QR) with ||res||=",
//                          rep.residual_norm);
        } else {
            // Identify suspect rows (largest residual entries).
            struct Pair { Precision val; Index row; };
            std::vector<Pair> rows; rows.reserve(set.m);
            for (Index i = 0; i < set.m; ++i) rows.push_back({ std::abs(resid[i]), i });
            const int K = std::min<int>(opt.suspect_rows_k, rows.size());
            std::partial_sort(rows.begin(), rows.begin() + K, rows.end(),
                              [](const Pair& a, const Pair& b){ return a.val > b.val; });

            rep.suspect_row_ids.clear();
            rep.suspect_row_ids.reserve(K);
            for (int k = 0; k < K; ++k) {
                const Index orig = set.kept_row_ids.empty() ? rows[k].row
                                                            : set.kept_row_ids[rows[k].row];
                rep.suspect_row_ids.push_back(orig);
            }
            logging::warning(false, "[ConstraintBuilder] Infeasible constraints (dense QR): ||C u - d|| = ",
                             rep.residual_norm, " (tol ", feas_tol, "). Top ", K, " suspect rows reported.");
        }
    }

    // 6) Summary log.
    {
        std::ostringstream oss;
//        oss << "[ConstraintBuilder] QR (dense fallback) finished: C(" << set.m << "x" << set.n << "), "
//            << "rank=" << rep.rank << ", redundant_rows~" << rep.n_redundant_rows << ". "
//            << (rep.homogeneous ? "Homogeneous constraints."
//                                : (std::string("Inhomogeneous: ||res||=") + std::to_string(rep.residual_norm) +
//                                   (rep.feasible ? " (feasible)." : " (INFEASIBLE!)")));
        rep.log = oss.str();
        logging::info(true, rep.log);
    }

    return {M, rep};
}

/******************************************************************************
 * @brief Build map with defaults (delegates to full overload).
 ******************************************************************************/
std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt; // defaults from the header
    return build(set, opt);
}

/******************************************************************************
 * @brief Build map with options; chooses dense or sparse path, assembles map & report.
 *
 * Dense path is preferred for small m (cheap, robust pivoting). For larger
 * constraint sets, a sparse QR is used with the same rank logic as the dense
 * path. Both paths share the same T/X assembly and reporting semantics.
 ******************************************************************************/
std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set, const Options& opt)
{
    // Heuristic: use dense QR for small numbers of rows (robust & inexpensive).
    if (set.m <= 256) {
        return build_dense_qr(set, opt);
    }

    // --- Sparse path (large m) ------------------------------------------------
    ConstraintMap M;
    Report rep;
    rep.m = set.m;
    rep.n = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    QR qr;
    qr.setPivotThreshold(0.0); // decide rank explicitly via our tolerance
    qr.compute(set.C);
    logging::error(qr.info() == Eigen::Success,
                   "[ConstraintBuilder] SparseQR failed to factorize C (", set.m, "x", set.n, ").");

    // Retrieve R and compress for predictable iterators.
    SparseMatrix R = qr.matrixR();
    R.makeCompressed();

    // Robust rank: scan ALL diagonals of the triangular block and count above tol.
    const int rmax = std::min<int>(set.m, set.n);
    Precision max_abs_diag = 0;
    std::vector<Precision> diag_abs(rmax, 0);
    for (int k = 0; k < rmax; ++k) {
        diag_abs[k] = read_diag_abs(R, k);
        if (diag_abs[k] > max_abs_diag) max_abs_diag = diag_abs[k];
    }
    rep.R11_max_diag = max_abs_diag;
    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);

    int r = 0;
    for (int k = 0; k < rmax; ++k) if (diag_abs[k] > tol) ++r; // no early break
    rep.rank = r;

    // Column permutation and assembly of (X, T).
    Eigen::VectorXi p = qr.colsPermutation().indices();
    SparseMatrix Rtop = R.topRows(r);
    build_T_from_R_and_perm(set.n, r, p, Rtop, M.slaves_, M.masters_, M.X_, M.T_);

    // Populate sizes & defaults; assess feasibility if d ≠ 0.
    M.n_   = set.n;
    M.r_   = r;
    M.nm_  = set.n - r;
    M.u_p_ = DynamicVector::Zero(set.n);
    rep.feasible = true;
    rep.d_norm   = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

    if (!rep.homogeneous) {
        // Least-squares particular solution (sparse QR).
        DynamicVector u_any = qr.solve(set.d);
        DynamicVector resid = set.C * u_any - set.d;
        rep.residual_norm   = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any;
            DynamicVector s  = DynamicVector::Zero(M.r_);
            for (int i = 0; i < M.r_; ++i) s[i] = u_any[M.slaves_[i]] - Xq[i];
            for (int i = 0; i < M.r_; ++i) M.u_p_[M.slaves_[i]] = s[i];

            logging::info(true, "[ConstraintBuilder] Built particular solution u_p (sparse QR) with ||res||=",
                          rep.residual_norm);
        } else {
            struct Pair { Precision val; Index row; };
            std::vector<Pair> rows; rows.reserve(set.m);
            for (Index i = 0; i < set.m; ++i) rows.push_back({ std::abs(resid[i]), i });
            const int K = std::min<int>(opt.suspect_rows_k, rows.size());
            std::partial_sort(rows.begin(), rows.begin() + K, rows.end(),
                              [](const Pair& a, const Pair& b){ return a.val > b.val; });

            rep.suspect_row_ids.clear();
            rep.suspect_row_ids.reserve(K);
            for (int k = 0; k < K; ++k) {
                const Index orig = set.kept_row_ids.empty() ? rows[k].row
                                                            : set.kept_row_ids[rows[k].row];
                rep.suspect_row_ids.push_back(orig);
            }
            logging::warning(false, "[ConstraintBuilder] Infeasible constraints (sparse QR): ||C u - d|| = ",
                             rep.residual_norm, " (tol ", feas_tol, "). Top ", K, " suspect rows reported.");
        }
    }

    // Summary
    {
        std::ostringstream oss;
        oss << "[ConstraintBuilder] QR (sparse) finished: C(" << set.m << "x" << set.n << "), "
            << "rank=" << rep.rank << ", redundant_rows~" << rep.n_redundant_rows << ". "
            << (rep.homogeneous ? "Homogeneous constraints."
                                : (std::string("Inhomogeneous: ||res||=") + std::to_string(rep.residual_norm) +
                                   (rep.feasible ? " (feasible)." : " (INFEASIBLE!)")));
        rep.log = oss.str();
        logging::info(true, rep.log);
    }

    return {M, rep};
}

} // namespace fem::constraint
