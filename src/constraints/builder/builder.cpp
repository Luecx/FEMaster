/******************************************************************************
 * @file builder.cpp
 * @brief Orchestrates modular constraint building pipeline.
 ******************************************************************************/

#include "builder.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../constraint_map.h"    // for ConstraintMap members used below
#include "./assemble_TX.h"
#include "./build_x.h"
#include "./detect_rank.h"
#include "./factorize_qr.h"
#include "./invariants.h"
#include "./particular_solution.h"
#include "./partition.h"
#include "./preprocess.h"
#include "./timings.h"

namespace fem::constraint {

std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt;
    return build(set, opt);
}

std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set,
                                                                             const Options&       opt) {
    using namespace std;

    ConstraintMap M;
    Report        rep {};
    rep.m = set.m;
    rep.n = set.n;

    // Timers
    Time t_pre = 0;
    Time t_qr = 0;
    Time t_rank = 0;
    Time t_x = 0;
    Time t_particular = 0;
    Time t_invar = 0;
    Time t_assemble = 0;
    Time t_partition = 0;

    // ---------------- 1) Preprocess ----------------
    PreprocessInput  pin {set.C,
                         set.d,
                         static_cast<int>(set.n),
                         static_cast<int>(set.m),
                         (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0))};
    PreprocessOutput pre;
    t_pre           = Timer::measure_time([&] { pre = preprocess_constraints(pin); });

    rep.homogeneous = pin.homogeneous;
    rep.d_norm      = (set.d.size() == 0) ? Precision(0) : set.d.norm();

    // -------- Early out: no used columns left → identity on non-fixed --------
    if (pre.n_use == 0) {
        std::vector<Index> masters_glob;
        masters_glob.reserve(set.n);
        for (Index j = 0; j < set.n; ++j)
            if (!pre.is_fixed_col[static_cast<int>(j)])
                masters_glob.push_back(j);

        M.n_       = set.n;
        M.r_       = 0;
        M.nm_      = (Index) masters_glob.size();
        M.masters_ = std::move(masters_glob);
        M.slaves_.clear();

        // T = identity on masters
        M.T_.resize(set.n, (int) M.nm_);
        {
            std::vector<Eigen::Triplet<Precision>> trips;
            trips.reserve((size_t) M.nm_);
            for (int j = 0; j < (int) M.nm_; ++j)
                trips.emplace_back(M.masters_[j], j, Precision(1));
            M.T_.setFromTriplets(trips.begin(), trips.end());
            M.T_.makeCompressed();
        }
        // X empty
        M.X_.resize(0, (int) M.nm_);

        // u_p: only fixed DOFs
        M.u_p_ = DynamicVector::Zero(set.n);
        for (int j = 0; j < (int) set.n; ++j)
            if (pre.is_fixed_col[j])
                M.u_p_[j] = pre.fixed_val[j];

        // Report
        int kept_rows = 0;
        for (char k : pre.keep_row)
            kept_rows += (k ? 1 : 0);
        rep.rank             = 0;
        rep.n_redundant_rows = kept_rows;
        rep.feasible         = true;
        rep.R11_max_diag     = 0;
        rep.residual_norm    = 0;

        // Timing table
        print_constraint_builder_timing({
            {"Preprocess", t_pre},
            {"QR", 0},
            {"Rank", 0},
            {"Build X", 0},
            {"Partition", 0},
            {"Assemble", 0},
#ifndef NDEBUG
            {"Invariants", 0},
#endif
            {"Particular", 0},
        });

        return {M, rep};
    }

    // ---------------- 2) QR ----------------
    QRSettings qrs {std::min<Precision>(0.1, std::max<Precision>(0, opt.rank_tol_rel))};
    QRResult   qr;
    // CHANGED: fill by reference to avoid moving Eigen::SparseQR
    t_qr = Timer::measure_time([&] { factorize_sparse_qr(pre.C_use, qrs, qr); });

    // ---------------- 3) Rank ----------------
    RankSettings rs {opt.rank_tol_rel};
    RankInfo     rk;
    t_rank           = Timer::measure_time([&] { rk = detect_rank_from_R(qr.R, rs); });
    rep.rank         = rk.r;
    rep.R11_max_diag = rk.max_abs_diag;

    // ---------------- 4) X columns ----------------
    XCols Xc;
    t_x              = Timer::measure_time([&] { Xc = build_X_cols_from_R(qr.R, rk.r); });
    const int nm_use = (int) pre.n_use - rk.r;
    (void) nm_use;

    // ---------------- 5) Partition + global mapping ----------------
    std::vector<char> is_used(set.n, 0);
    for (int c : pre.used)
        is_used[c] = 1;
    PartitionInput  part_in {qr.P, rk.r, &pre.used, &pre.is_fixed_col, (int) set.n};
    PartitionOutput part;
    t_partition = Timer::measure_time([&] { part = partition_and_map(part_in, is_used); });

    // ---------------- 6) Assemble T and X ----------------
    AssembleInput  ain {(int) set.n, rk.r, part.slaves_loc, part.masters_loc, pre.used, part.masters_glob};
    AssembleOutput ax;
    t_assemble = Timer::measure_time([&] { ax = assemble_T_and_X(ain, Xc); });

    // Populate map
    M.n_       = set.n;
    M.r_       = rk.r;
    M.nm_      = (Index) part.masters_glob.size();
    M.masters_ = std::move(part.masters_glob);
    M.slaves_  = std::move(part.slaves_glob);
    M.T_       = std::move(ax.T);
    M.X_       = std::move(ax.X);
    M.u_p_     = DynamicVector::Zero(set.n);
    for (int j = 0; j < (int) set.n; ++j)
        if (pre.is_fixed_col[j])
            M.u_p_[j] = pre.fixed_val[j];

#ifndef NDEBUG
    t_invar = Timer::measure_time([&] { check_invariants(M, set.C); });
#endif

    // ---------------- 7) Particular solution (inhom) ----------------
    t_particular = Timer::measure_time([&] {
        if (!rep.homogeneous) {
            ParticularInput
                 pin2 {rep.homogeneous, pre.C_use, pre.d_mod, pre.used, pre.is_fixed_col, pre.fixed_val, set.d};
            auto pout         = compute_particular_and_project(pin2, qr.qr, M, opt.feas_tol_rel, rep.d_norm);
            M.u_p_            = std::move(pout.u_p);
            rep.residual_norm = pout.residual_norm;
            rep.feasible      = pout.feasible;
        } else {
            rep.residual_norm = 0;
            rep.feasible      = true;
        }
    });

    // Redundant rows
    int kept_rows = 0;
    for (char k : pre.keep_row)
        kept_rows += (k ? 1 : 0);
    rep.n_redundant_rows = (kept_rows >= rep.rank) ? (kept_rows - rep.rank) : 0;

    // Timing table
    print_constraint_builder_timing({
        {"Preprocess", t_pre},
        {"QR", t_qr},
        {"Rank", t_rank},
        {"Build X", t_x},
        {"Partition", t_partition},
        {"Assemble", t_assemble},
#ifndef NDEBUG
        {"Invariants", t_invar},
#endif
        {"Particular", t_particular},
    });

    return {M, rep};
}

}    // namespace fem::constraint
