
#include "constraint_builder.h"
#include "../constraint_map.h"
#include "../../core/logging.h"
#include "../../core/timer.h"

// step headers
#include "./cb_memory.h"
#include "./cb_types.h"
#include "./cb_scan_simple.h"
#include "./cb_substitute.h"
#include "./cb_compress.h"
#include "./cb_qr_rank.h"
#include "./cb_build_x_sparse.h"
#include "./cb_partition.h"
#include "./cb_build_TX.h"
#include "./cb_particular.h"

namespace fem { namespace constraint {

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt; return build(set, opt);
}

std::pair<ConstraintMap, ConstraintBuilder::Report>
ConstraintBuilder::build(const ConstraintSet& set, const Options& opt)
{
    using Perm = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;

    struct Times { Time scan_simple{0}, apply_subst{0}, compress{0}, qr{0}, rank_detect{0},
                          build_Xd{0}, partition{0}, build_T{0}, build_X{0},
                          invariants{0}, particular{0}, summary_log{0}; } t;

    MemoryTracker mem; // live objects accounting

    ConstraintMap M;
    Report rep; rep.m=set.m; rep.n=set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));
    rep.d_norm = (set.d.size() == 0) ? Precision(0) : set.d.norm();

    // 1) Scan simple rows
    std::vector<SimpleRow> singles;
    t.scan_simple = Timer::measure_time([&]{ singles = scan_single_nnz_rows(set.C, set.d, rep.homogeneous); });

    // 2) Substitute fixed
    SubstituteResult sub;
    t.apply_subst = Timer::measure_time([&]{ sub = substitute_fixed(set.C, set.d, singles); });

    // Track memory of updated C and d
    mem.add_sparse("C_after_substitute", sub.C);
    mem.add_vector("d_after_substitute", std::vector<Precision>((size_t)sub.d.size()));

    // 3) Compress zero columns with row filter
    std::vector<int> used; Eigen::VectorXi old2new; SparseMatrix C_use;
    t.compress = Timer::measure_time([&]{ compress_zero_columns_with_row_filter(sub.C, sub.keep_row, used, old2new, C_use); });
    mem.add_sparse("C_use", C_use);

    if ((int)used.size() == 0) {
        // identity map (except fixed cols in u_p)
        std::vector<Index> masters_glob; masters_glob.reserve(set.n);
        for (int j = 0; j < set.n; ++j) if (!sub.is_fixed_col[j]) masters_glob.push_back(j);
        M.n_ = set.n; M.r_ = 0; M.nm_ = (Index)masters_glob.size(); M.masters_ = std::move(masters_glob);
        M.slaves_.clear();
        M.T_ = SparseMatrix(set.n, (int)M.nm_); M.T_.setIdentity();
        M.X_.resize(0, (int)M.nm_);
        M.u_p_ = DynamicVector::Zero(set.n);
        for (int j = 0; j < set.n; ++j) if (sub.is_fixed_col[j]) M.u_p_[j] = sub.fixed_val[j];

        // memory accounting for outputs
        mem.add_sparse("T_full", M.T_);
        mem.add_sparse("X", M.X_);
        mem.add_vector("u_p", std::vector<Precision>((size_t)M.u_p_.size()));
        mem.add_vector("masters", M.masters_);
        mem.add_vector("slaves", M.slaves_);

        // timings + memory
        t.summary_log = Timer::measure_time([&]{
            logging::info(true, "");
            logging::info(true, "ConstraintBuilder timings (ms):");
            logging::info(true, "  Scan simple rows:          ", std::setw(8), t.scan_simple);
            logging::info(true, "  Apply substitutions:       ", std::setw(8), t.apply_subst);
            logging::info(true, "  Compress columns:          ", std::setw(8), t.compress);
            logging::info(true, "  (No QR path taken)");
            MemoryTracker::log_summary(mem);
        });
        rep.rank = 0; rep.n_redundant_rows = 0; rep.feasible = true; rep.residual_norm = 0;
        return {M, rep};
    }

    // 4/5) QR + rank
    QRRank qr;
    t.qr = Timer::measure_time([&]{
        factor_and_rank_inplace(C_use, opt.rank_tol_rel, qr);
    });
    logging::error(qr.sqr.info() == Eigen::Success, "[ConstraintBuilder] SparseQR failed.");
    mem.add_sparse("R", qr.R);
    rep.R11_max_diag = qr.R11_max_diag;
    rep.rank = qr.r;


    const int r = qr.r; const int nm_use = (int)used.size() - r;

    // 6) Build X (sparse back-sub)
    XBuildOut xb;
    t.build_Xd = Timer::measure_time([&]{ xb = build_X_sparse(qr.R, r, nm_use); });
    mem.add_nested_vector_pair("X_cols", xb.X_cols);
    mem.add_nested_vector("R11_rows", xb.R11_rows);
    mem.add_vector("R11_diag", xb.R11_diag);

    // 7) Partition by permutation
    PartitionOut part;
    t.partition = Timer::measure_time([&]{ part = partition_by_P(qr.sqr.colsPermutation(), r, nm_use); });

    // mark used columns
    std::vector<char> is_used(set.n, 0); for (int c : used) is_used[c] = 1;

    // 8/9) Build T and X
    TXOut tx;
    t.build_T = Timer::measure_time([&]{
        tx = build_T_and_X(set.n, r, used, part.slaves_loc, part.masters_loc, is_used,
                           sub.is_fixed_col, xb.X_cols);
    });
    t.build_X = 0; // already built inside

    mem.add_sparse("T_full", tx.T_full);
    mem.add_sparse("X", tx.X);
    mem.add_vector("masters", tx.masters_glob);
    mem.add_vector("slaves", tx.slaves_glob);

    // populate map
    M.n_ = set.n; M.r_ = r; M.nm_ = (Index)tx.masters_glob.size();
    M.masters_ = std::move(tx.masters_glob);
    M.slaves_  = std::move(tx.slaves_glob);
    M.T_       = std::move(tx.T_full);
    M.X_       = std::move(tx.X);
    M.u_p_     = DynamicVector::Zero(set.n);
    for (int j = 0; j < set.n; ++j) if (sub.is_fixed_col[j]) M.u_p_[j] = sub.fixed_val[j];

#ifndef NDEBUG
    // Invariants
    t.invariants = Timer::measure_time([&]{
        const int nm = (int)M.nm_;
        DynamicVector e = DynamicVector::Zero(nm);
        double max_norm = 0.0;
        for (int k = 0; k < nm; ++k) {
            e.setZero(); e[k] = 1;
            DynamicVector cu = set.C * (M.T_ * e);
            max_norm = std::max(max_norm, (double)cu.norm());
        }
        logging::error(max_norm <= 1e-12, "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ", max_norm);
        DynamicVector q = DynamicVector::Random(M.nm_);
        DynamicVector u_explicit = M.T_ * q + M.u_p_;
        DynamicVector u_fast     = DynamicVector::Zero(M.n_);
        M.apply_T(q, u_fast); u_fast.array() += M.u_p_.array();
        logging::error( (u_explicit - u_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_T mismatch explicit T.");
        DynamicVector y = DynamicVector::Random(M.n_);
        DynamicVector z_explicit = M.T_.transpose() * y;
        DynamicVector z_fast     = DynamicVector::Zero(M.nm_);
        M.apply_Tt(y, z_fast);
        logging::error( (z_explicit - z_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    });
#endif

    // 10) Particular solution (if inhomogeneous)
    ParticularOut p;
    t.particular = Timer::measure_time([&]{
        p = compute_particular(C_use, qr.sqr, set.C, set.d, used,
                               sub.is_fixed_col, sub.fixed_val, set.n,
                               opt.feas_tol_rel, rep.d_norm, M.X_, M.masters_, M.slaves_);
    });

    rep.feasible = p.feasible; rep.residual_norm = p.residual_norm;
    if (p.feasible) M.u_p_ = std::move(p.u_p);

    mem.add_vector("u_p", std::vector<Precision>((size_t)M.u_p_.size()));

    // Redundant rows estimate
    int kept_rows = 0; for (char k : sub.keep_row) kept_rows += (k ? 1 : 0);
    rep.n_redundant_rows = (kept_rows >= rep.rank) ? (kept_rows - rep.rank) : 0;

    // Summary timings + memory
    t.summary_log = Timer::measure_time([&]{
        logging::info(true, "");
        logging::info(true, "ConstraintBuilder timings (ms):");
        logging::info(true, "  Scan simple rows:          ", std::setw(8), t.scan_simple);
        logging::info(true, "  Apply substitutions:       ", std::setw(8), t.apply_subst);
        logging::info(true, "  Compress columns:          ", std::setw(8), t.compress);
        logging::info(true, "  Solve QR:                  ", std::setw(8), t.qr);
        logging::info(true, "  Rank detection:            ", std::setw(8), t.rank_detect);
        logging::info(true, "  Build X (sparse solve):    ", std::setw(8), t.build_Xd);
        logging::info(true, "  Partition permutation:     ", std::setw(8), t.partition);
        logging::info(true, "  Build T/X (sparse):        ", std::setw(8), t.build_T);
#ifndef NDEBUG
        logging::info(true, "  Invariants check:          ", std::setw(8), t.invariants);
#endif
        logging::info(true, "  Particular solution:       ", std::setw(8), t.particular);
        // MemoryTracker::log_summary(mem);
    });

    return {M, rep};
}

}} // namespace fem::constraint