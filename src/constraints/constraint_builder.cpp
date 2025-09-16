/******************************************************************************
 * @file constraint_builder.cpp
 * @brief Implementation: build a reduced null-space map u = u_p + T q from C u = d.
 ******************************************************************************/

#include "constraint_builder.h"
#include "constraint_map.h"
#include "../core/logging.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <Eigen/QR>   // ColPivHouseholderQR
#include <algorithm>
#include <sstream>
#include <unordered_map>

namespace fem::constraint {

using QR = Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>;

/*** Utilities ***************************************************************/

// Read |R(k,k)| from sparse upper-triangular column k
static Precision read_diag_abs(const SparseMatrix& R, int k) {
    for (SparseMatrix::InnerIterator it(R, k); it; ++it) {
        if (it.row() == k) return std::abs(it.value());
        if (it.row() >  k) break;
    }
    return Precision(0);
}

// Build list of used (nonzero) columns and a compressed C_use = C(:, used)
static void compress_zero_columns(const SparseMatrix& C,
                                  std::vector<int>&   used_cols,
                                  Eigen::VectorXi&    old2new,
                                  SparseMatrix&       C_use)
{
    const int n = C.cols();
    used_cols.clear();
    used_cols.reserve(n);

    // Identify columns with at least one nonzero
    for (int j = 0; j < C.outerSize(); ++j) {
        bool nonzero = false;
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) { nonzero = true; break; }
        if (nonzero) used_cols.push_back(j);
    }

    // Map: old col -> new col (or -1 if unused)
    old2new = Eigen::VectorXi::Constant(n, -1);
    for (int k = 0; k < (int)used_cols.size(); ++k) old2new[used_cols[k]] = k;

    // Build C_use by copying only used columns
    C_use.resize(C.rows(), (int)used_cols.size());
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C.nonZeros());
    for (int jnew = 0; jnew < (int)used_cols.size(); ++jnew) {
        int j = used_cols[jnew];
        for (SparseMatrix::InnerIterator it(C, j); it; ++it) {
            trips.emplace_back(it.row(), jnew, it.value());
        }
    }
    C_use.setFromTriplets(trips.begin(), trips.end());
    C_use.makeCompressed();
}

// Build X (r x nm_use) and local partitions from R and permutation (local n_use)
static void build_X_and_partition_from_R(
    int n_use, int r, const Eigen::VectorXi& permCols, const SparseMatrix& Rtop,
    std::vector<int>& slaves_loc, std::vector<int>& masters_loc, SparseMatrix& X_use)
{
    const int nm_use = n_use - r;

    // Partition columns in local (compressed) index space
    slaves_loc.resize(r);
    masters_loc.resize(nm_use);
    for (int k = 0; k < r; ++k) slaves_loc[k] = permCols[k];
    for (int k = r; k < n_use; ++k) masters_loc[k - r] = permCols[k];

    // Extract R11 (r x r) and R12 (r x nm_use)
    SparseMatrix R11 = Rtop.leftCols(r);
    SparseMatrix R12 = Rtop.rightCols(nm_use);

    // Solve X_use = - R11^{-1} R12 (upper-triangular)
    X_use.resize(r, nm_use);
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
    X_use.setFromTriplets(xTrips.begin(), xTrips.end());
    X_use.makeCompressed();
}

// Build full-size T (n x nm_total) and global master/slave lists.
// masters_loc/slaves_loc are in compressed space; 'used' lifts to global DOFs.
static void assemble_global_T_from_local(
    int n_full,
    const std::vector<int>& used,
    const std::vector<int>& slaves_loc,
    const std::vector<int>& masters_loc,
    const SparseMatrix&     X_use,
    std::vector<Index>&     slaves_glob,
    std::vector<Index>&     masters_glob,
    SparseMatrix&           T_full)
{
    const int r = (int)slaves_loc.size();
    const int nm_use = (int)masters_loc.size();

    // Map local -> global DOFs
    slaves_glob.resize(r);
    for (int i = 0; i < r; ++i) slaves_glob[i] = used[slaves_loc[i]];

    // Start with masters that were in 'used'
    masters_glob.resize(0);
    masters_glob.reserve(n_full - r);
    for (int j = 0; j < nm_use; ++j) masters_glob.push_back(used[masters_loc[j]]);

    // Append all never-used columns as additional masters
    // Build a marker for used cols
    std::vector<char> is_used(n_full, 0);
    for (int c : used) is_used[c] = 1;

    // Also mark the already-added masters (from used) to avoid duplicates (not strictly needed)
    std::vector<char> is_master(n_full, 0);
    for (Index g : masters_glob) is_master[(int)g] = 1;

    for (int g = 0; g < n_full; ++g) {
        if (!is_used[g]) {
            masters_glob.push_back(g);
        }
    }

    const int nm_total = (int)masters_glob.size();

    // Map: global master DOF -> column index in T_full
    std::unordered_map<int,int> master_col_of_dof;
    master_col_of_dof.reserve(nm_total * 2);
    for (int j = 0; j < nm_total; ++j) master_col_of_dof[(int)masters_glob[j]] = j;

    // Assemble T_full: identity on masters; scatter X_use into slave rows
    T_full.resize(n_full, nm_total);
    std::vector<Eigen::Triplet<Precision>> tTrips;
    tTrips.reserve(nm_total + X_use.nonZeros());

    // Identity on masters
    for (int j = 0; j < nm_total; ++j) {
        tTrips.emplace_back((int)masters_glob[j], j, Precision(1));
    }

    // Scatter X_use (which only references masters from 'used')
    for (int jloc = 0; jloc < X_use.outerSize(); ++jloc) {
        // Global master dof corresponding to this local master column
        int g_master = used[masters_loc[jloc]];
        int jglob = master_col_of_dof[g_master];

        for (SparseMatrix::InnerIterator it(X_use, jloc); it; ++it) {
            int i = it.row();               // slave index in local list
            int g_slave = (int)slaves_glob[i];
            tTrips.emplace_back(g_slave, jglob, it.value());
        }
    }

    T_full.setFromTriplets(tTrips.begin(), tTrips.end());
    T_full.makeCompressed();
}

// Expand a reduced solution u_use (size n_use) to full size n_full using 'used'
static DynamicVector expand_solution_to_full(const DynamicVector& u_use,
                                             const std::vector<int>& used,
                                             int n_full)
{
    DynamicVector u_full = DynamicVector::Zero(n_full);
    for (int k = 0; k < (int)used.size(); ++k) {
        u_full[used[k]] = u_use[k];
    }
    return u_full;
}

/*** Dense path ***************************************************************/

static std::pair<ConstraintMap, ConstraintBuilder::Report>
build_dense_qr(const ConstraintSet& set, const ConstraintBuilder::Options& opt)
{
    ConstraintMap M;
    ConstraintBuilder::Report rep;

    rep.m = set.m;
    rep.n = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    // Compress zero columns
    std::vector<int> used;
    Eigen::VectorXi old2new; // not used on dense path, but kept for symmetry
    SparseMatrix C_use;
    compress_zero_columns(set.C, used, old2new, C_use);

    const int n_use = (int)used.size();

    // If nothing is used, everything is a master; trivial map.
    if (n_use == 0) {
        M.n_  = set.n;
        M.r_  = 0;
        M.nm_ = set.n;
        M.masters_.resize(set.n);
        M.slaves_.clear();
        for (int j = 0; j < set.n; ++j) M.masters_[j] = j;

        // Build T = identity
        M.T_ = SparseMatrix(set.n, set.n);
        M.T_.setIdentity();
        M.u_p_ = DynamicVector::Zero(set.n);

        rep.rank = 0;
        rep.n_redundant_rows = rep.m; // all rows redundant if any
        rep.feasible = true;
        rep.d_norm = (set.d.size() == 0) ? 0 : set.d.norm();
        return {M, rep};
    }

    // Dense QR on C_use
    Eigen::MatrixXd Cd = Eigen::MatrixXd(C_use);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dqr(Cd);

    Eigen::MatrixXd R = dqr.matrixR()
                            .topLeftCorner(set.m, n_use)
                            .template triangularView<Eigen::Upper>();

    // Rank decision
    Precision max_abs_diag = 0;
    const int rmax = std::min<int>(set.m, n_use);
    for (int k = 0; k < rmax; ++k) max_abs_diag = std::max<Precision>(max_abs_diag, std::abs(R(k,k)));
    rep.R11_max_diag = max_abs_diag;
    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);

    int r = 0;
    for (int k = 0; k < rmax; ++k) if (std::abs(R(k,k)) > tol) ++r;
    rep.rank = r;

    // Convert top r rows of R to sparse Rs (r x n_use)
    SparseMatrix Rs(r, n_use);
    std::vector<Eigen::Triplet<Precision>> trips;
    trips.reserve(C_use.nonZeros());
    for (int i = 0; i < r; ++i) {
        for (int j = i; j < n_use; ++j) {
            double v = R(i, j);
            if (v != 0.0) trips.emplace_back(i, j, v);
        }
    }
    Rs.setFromTriplets(trips.begin(), trips.end());
    Rs.makeCompressed();

    // Build X (local) and partitions (local)
    std::vector<int> slaves_loc, masters_loc;
    SparseMatrix X_use;
    build_X_and_partition_from_R(n_use, r, dqr.colsPermutation().indices(), Rs,
                                 slaves_loc, masters_loc, X_use);

    // Lift to global, assemble full T
    assemble_global_T_from_local(set.n, used, slaves_loc, masters_loc,
                                 X_use, M.slaves_, M.masters_, M.T_);

    // Populate sizes and defaults
    M.n_   = set.n;
    M.r_   = r;
    M.nm_  = (Index)M.masters_.size();
    M.u_p_ = DynamicVector::Zero(set.n);

    rep.feasible = true;
    rep.d_norm   = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

    if (!rep.homogeneous) {
        // Least-squares solution on C_use, then expand
        Eigen::VectorXd u_use = dqr.solve(set.d);
        DynamicVector u_any = expand_solution_to_full(u_use, used, set.n);

        DynamicVector resid = set.C * u_any - set.d;
        rep.residual_norm   = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            // Build u_p as in the original code
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any; // M.X_ not yet set — set it now from X_use

            // We need M.X_ to match slave rows and full master columns order.
            // Build M.X_ as r x nm_total by scattering X_use into the correct master columns.
            {
                M.X_.resize((int)M.slaves_.size(), (int)M.masters_.size());
                std::vector<Eigen::Triplet<Precision>> xTrips2;
                xTrips2.reserve(X_use.nonZeros());
                // Map global master dof -> column index
                std::unordered_map<int,int> master_col_of_dof;
                for (int j = 0; j < (int)M.masters_.size(); ++j)
                    master_col_of_dof[(int)M.masters_[j]] = j;

                for (int jloc = 0; jloc < X_use.outerSize(); ++jloc) {
                    int g_master = used[masters_loc[jloc]];
                    int jglob = master_col_of_dof[g_master];
                    for (SparseMatrix::InnerIterator it(X_use, jloc); it; ++it) {
                        int i = it.row();
                        xTrips2.emplace_back(i, jglob, it.value());
                    }
                }
                M.X_.setFromTriplets(xTrips2.begin(), xTrips2.end());
                M.X_.makeCompressed();
            }

            DynamicVector Xq2 = M.X_ * q_any;
            DynamicVector s  = DynamicVector::Zero(M.r_);
            for (int i = 0; i < M.r_; ++i) s[i] = u_any[M.slaves_[i]] - Xq2[i];
            for (int i = 0; i < M.r_; ++i) M.u_p_[M.slaves_[i]] = s[i];
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
            logging::warning(false, "[ConstraintBuilder] Infeasible constraints (dense QR): ||C u - d|| = ",
                             rep.residual_norm, " (tol ", feas_tol, "). Top ", K, " suspect rows reported.");
        }
    } else {
        // Build M.X_ consistent with M.T_ (X only on masters from 'used'; zeros elsewhere)
        M.X_.resize((int)M.slaves_.size(), (int)M.masters_.size());
        std::vector<Eigen::Triplet<Precision>> xTrips2;
        // Map global master dof -> column index
        std::unordered_map<int,int> master_col_of_dof;
        for (int j = 0; j < (int)M.masters_.size(); ++j)
            master_col_of_dof[(int)M.masters_[j]] = j;

        for (int jloc = 0; jloc < X_use.outerSize(); ++jloc) {
            int g_master = used[masters_loc[jloc]];
            int jglob = master_col_of_dof[g_master];
            for (SparseMatrix::InnerIterator it(X_use, jloc); it; ++it) {
                int i = it.row();
                xTrips2.emplace_back(i, jglob, it.value());
            }
        }
        M.X_.setFromTriplets(xTrips2.begin(), xTrips2.end());
        M.X_.makeCompressed();
    }

    return {M, rep};
}

/*** Sparse path **************************************************************/

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt;
    return build(set, opt);
}

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set, const Options& opt)
{
    // Dense fallback for small m
    if (set.m <= 256) {
        return build_dense_qr(set, opt);
    }

    // Compress zero columns
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

    // If nothing is used, everything is a master; trivial map.
    if (n_use == 0) {
        M.n_  = set.n;
        M.r_  = 0;
        M.nm_ = set.n;
        M.masters_.resize(set.n);
        M.slaves_.clear();
        for (int j = 0; j < set.n; ++j) M.masters_[j] = j;

        M.T_ = SparseMatrix(set.n, set.n);
        M.T_.setIdentity();
        M.X_.resize(0, set.n);
        M.u_p_ = DynamicVector::Zero(set.n);

        rep.rank = 0;
        rep.n_redundant_rows = rep.m;
        rep.feasible = true;
        rep.d_norm = (set.d.size() == 0) ? 0 : set.d.norm();
        return {M, rep};
    }

    // Sparse QR on C_use
    QR qr;
    qr.setPivotThreshold(0.0);
    qr.compute(C_use);
    logging::error(qr.info() == Eigen::Success,
                   "[ConstraintBuilder] SparseQR failed to factorize C_use (",
                   set.m, "x", n_use, ").");

    // Retrieve R and analyze diag
    SparseMatrix R = qr.matrixR();
    R.makeCompressed();

    const int rmax = std::min<int>(set.m, n_use);
    Precision max_abs_diag = 0;
    std::vector<Precision> diag_abs(rmax, 0);
    for (int k = 0; k < rmax; ++k) {
        diag_abs[k] = read_diag_abs(R, k);
        if (diag_abs[k] > max_abs_diag) max_abs_diag = diag_abs[k];
    }
    rep.R11_max_diag = max_abs_diag;
    const Precision tol = std::max(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);

    int r = 0;
    for (int k = 0; k < rmax; ++k) if (diag_abs[k] > tol) ++r;
    rep.rank = r;

    // Build X in local space
    Eigen::VectorXi p = qr.colsPermutation().indices();
    SparseMatrix Rtop = R.topRows(r);
    std::vector<int> slaves_loc, masters_loc;
    SparseMatrix X_use;
    build_X_and_partition_from_R(n_use, r, p, Rtop, slaves_loc, masters_loc, X_use);

    // Lift to global, assemble full T
    assemble_global_T_from_local(set.n, used, slaves_loc, masters_loc,
                                 X_use, M.slaves_, M.masters_, M.T_);

    // Populate sizes
    M.n_   = set.n;
    M.r_   = r;
    M.nm_  = (Index)M.masters_.size();
    M.u_p_ = DynamicVector::Zero(set.n);

    rep.feasible = true;
    rep.d_norm   = (set.d.size() == 0) ? 0 : set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank) ? (rep.m - rep.rank) : 0;

    // Build M.X_ aligned to global master ordering (scatter X_use columns)
    {
        M.X_.resize((int)M.slaves_.size(), (int)M.masters_.size());
        std::vector<Eigen::Triplet<Precision>> xTrips2;
        xTrips2.reserve(X_use.nonZeros());
        std::unordered_map<int,int> master_col_of_dof;
        for (int j = 0; j < (int)M.masters_.size(); ++j)
            master_col_of_dof[(int)M.masters_[j]] = j;

        for (int jloc = 0; jloc < X_use.outerSize(); ++jloc) {
            int g_master = used[masters_loc[jloc]];
            int jglob = master_col_of_dof[g_master];
            for (SparseMatrix::InnerIterator it(X_use, jloc); it; ++it) {
                int i = it.row();
                xTrips2.emplace_back(i, jglob, it.value());
            }
        }
        M.X_.setFromTriplets(xTrips2.begin(), xTrips2.end());
        M.X_.makeCompressed();
    }

    if (!rep.homogeneous) {
        // Solve least-squares on compressed system, expand, and compute u_p
        DynamicVector u_use = qr.solve(set.d);
        DynamicVector u_any = expand_solution_to_full(u_use, used, set.n);

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

    return {M, rep};
}

} // namespace fem::constraint
