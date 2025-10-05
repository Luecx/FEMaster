/******************************************************************************
 * @file builder.cpp
 * @brief Implements the constraint builder pipeline.
 *
 * The builder orchestrates preprocessing, sparse QR factorisation, rank
 * detection, and assembly of the null-space transformation matrices.
 *
 * @see src/constraints/builder/builder.h
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "builder.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../constraint_map.h"
#include "assemble_TX.h"
#include "build_x.h"
#include "detect_rank.h"
#include "factorize_qr.h"
#include "invariants.h"
#include "particular_solution.h"
#include "partition.h"
#include "preprocess.h"
#include "timings.h"

namespace fem {
namespace constraint {

/******************************************************************************
 * @copydoc ConstraintBuilder::build(const ConstraintSet&)
 ******************************************************************************/
std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set) {
    Options opt;
    return build(set, opt);
}

/******************************************************************************
 * @copydoc ConstraintBuilder::build(const ConstraintSet&,const Options&)
 ******************************************************************************/
std::pair<ConstraintMap, ConstraintBuilder::Report> ConstraintBuilder::build(const ConstraintSet& set,
                                                                             const Options& opt) {
    ConstraintMap map;
    Report report;
    report.m = set.m;
    report.n = set.n;

    Time t_pre = 0;
    Time t_qr = 0;
    Time t_rank = 0;
    Time t_x = 0;
    Time t_particular = 0;
    Time t_invar = 0;
    Time t_assemble = 0;
    Time t_partition = 0;

    PreprocessInput pin{set.C,
                        set.d,
                        static_cast<int>(set.n),
                        static_cast<int>(set.m),
                        (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0))};
    PreprocessOutput pre;
    t_pre = Timer::measure_time([&] { pre = preprocess_constraints(pin); });

    report.homogeneous = pin.homogeneous;
    report.d_norm = (set.d.size() == 0) ? Precision(0) : set.d.norm();

    if (pre.n_use == 0) {
        std::vector<Index> masters_glob;
        masters_glob.reserve(set.n);
        for (Index j = 0; j < set.n; ++j) {
            if (!pre.is_fixed_col[static_cast<int>(j)]) {
                masters_glob.push_back(j);
            }
        }

        map.n_ = set.n;
        map.r_ = 0;
        map.nm_ = static_cast<Index>(masters_glob.size());
        map.masters_ = std::move(masters_glob);
        map.slaves_.clear();

        map.T_.resize(set.n, static_cast<int>(map.nm_));
        {
            std::vector<Eigen::Triplet<Precision>> trips;
            trips.reserve(static_cast<std::size_t>(map.nm_));
            for (int j = 0; j < static_cast<int>(map.nm_); ++j) {
                trips.emplace_back(map.masters_[j], j, Precision(1));
            }
            map.T_.setFromTriplets(trips.begin(), trips.end());
            map.T_.makeCompressed();
        }
        map.X_.resize(0, static_cast<int>(map.nm_));

        map.u_p_ = DynamicVector::Zero(set.n);
        for (int j = 0; j < static_cast<int>(set.n); ++j) {
            if (pre.is_fixed_col[j]) {
                map.u_p_[j] = pre.fixed_val[j];
            }
        }

        int kept_rows = 0;
        for (char keep : pre.keep_row) {
            kept_rows += (keep ? 1 : 0);
        }
        report.rank = 0;
        report.n_redundant_rows = kept_rows;
        report.feasible = true;
        report.R11_max_diag = 0;
        report.residual_norm = 0;

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

        return {map, report};
    }

    QRSettings qrs{std::min<Precision>(0.1, std::max<Precision>(0, opt.rank_tol_rel))};
    QRResult qr;
    t_qr = Timer::measure_time([&] { factorize_sparse_qr(pre.C_use, qrs, qr); });

    RankSettings rs{opt.rank_tol_rel};
    RankInfo rk;
    t_rank = Timer::measure_time([&] { rk = detect_rank_from_R(qr.R, rs); });
    report.rank = rk.r;
    report.R11_max_diag = rk.max_abs_diag;

    XCols xcols;
    t_x = Timer::measure_time([&] { xcols = build_X_cols_from_R(qr.R, rk.r); });

    std::vector<char> is_used(set.n, 0);
    for (int column : pre.used) {
        is_used[column] = 1;
    }
    PartitionInput part_in{qr.P, rk.r, &pre.used, &pre.is_fixed_col, static_cast<int>(set.n)};
    PartitionOutput part;
    t_partition = Timer::measure_time([&] { part = partition_and_map(part_in, is_used); });

    AssembleInput assemble_in{static_cast<int>(set.n),
                              rk.r,
                              part.slaves_loc,
                              part.masters_loc,
                              pre.used,
                              part.masters_glob};
    AssembleOutput assemble_out;
    t_assemble = Timer::measure_time([&] { assemble_out = assemble_T_and_X(assemble_in, xcols); });

    map.n_ = set.n;
    map.r_ = rk.r;
    map.nm_ = static_cast<Index>(part.masters_glob.size());
    map.masters_ = std::move(part.masters_glob);
    map.slaves_ = std::move(part.slaves_glob);
    map.T_ = std::move(assemble_out.T);
    map.X_ = std::move(assemble_out.X);
    map.u_p_ = DynamicVector::Zero(set.n);
    for (int j = 0; j < static_cast<int>(set.n); ++j) {
        if (pre.is_fixed_col[j]) {
            map.u_p_[j] = pre.fixed_val[j];
        }
    }

#ifndef NDEBUG
    t_invar = Timer::measure_time([&] { check_invariants(map, set.C); });
#endif

    t_particular = Timer::measure_time([&] {
        if (!report.homogeneous) {
            ParticularInput pin2{report.homogeneous,
                                 pre.C_use,
                                 pre.d_mod,
                                 pre.used,
                                 pre.is_fixed_col,
                                 pre.fixed_val,
                                 set.d};
            auto pout = compute_particular_and_project(pin2, qr.qr, map, opt.feas_tol_rel, report.d_norm);
            map.u_p_ = std::move(pout.u_p);
            report.residual_norm = pout.residual_norm;
            report.feasible = pout.feasible;
        } else {
            report.residual_norm = 0;
            report.feasible = true;
        }
    });

    int kept_rows = 0;
    for (char keep : pre.keep_row) {
        kept_rows += (keep ? 1 : 0);
    }
    report.n_redundant_rows = (kept_rows >= report.rank) ? (kept_rows - report.rank) : 0;

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

    return {map, report};
}

} // namespace constraint
} // namespace fem
