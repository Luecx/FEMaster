/******************************************************************************
 * @file constraint_builder.cpp
 * @brief Build a reduced null-space map u = u_p + T q from C u = d.
 *
 * SparseQR-only pipeline (COLAMD ordering, single-thread friendly):
 *   - Compress zero columns of C → C_use
 *   - SparseQR: C_use * P = Q * R   (pivotThreshold=0; rank via diag(R))
 *   - Rank r
 *   - Solve X = -R11^{-1} R12
 *   - Optional drop tolerance on X
 *   - Build T and X as sparse
 *   - u_p via SparseQR::solve if inhomogeneous
 *
 * @date    17.09.2025  (timed version, COLAMD fixed)
 * @author  Finn
 ******************************************************************************/

#include "constraint_builder.h"
#include "constraint_map.h"
#include "../core/logging.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <chrono>

namespace fem::constraint {

using Clock = std::chrono::high_resolution_clock;

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

    auto t_start = Clock::now();
    auto lap = [&](const char* msg, Clock::time_point& t0) {
        auto t1 = Clock::now();
        double ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        logging::info(true, msg, " : ", ms, " ms");
        t0 = t1;
    };

    // ---- Compress columns ---------------------------------------------------
    std::vector<int> used;
    Eigen::VectorXi old2new;
    SparseMatrix C_use;
    compress_zero_columns(set.C, used, old2new, C_use);
    lap("compress_zero_columns", t_start);

    ConstraintMap M;
    Report rep;
    rep.m = set.m;
    rep.n = set.n;
    rep.homogeneous = (set.d.size() == 0) || (set.d.lpNorm<Eigen::Infinity>() == Precision(0));

    const int n_use = (int)used.size();

    if (n_use == 0) {
        M.n_  = set.n; M.r_ = 0; M.nm_ = set.n;
        M.masters_.resize(set.n);
        for (int j = 0; j < set.n; ++j) M.masters_[j] = j;
        M.slaves_.clear();
        M.T_ = SparseMatrix(set.n, set.n); M.T_.setIdentity();
        M.X_.resize(0, set.n);
        M.u_p_ = DynamicVector::Zero(set.n);
        rep.rank = 0;
        rep.n_redundant_rows = rep.m;
        rep.feasible = true;
        rep.d_norm = (set.d.size() == 0) ? 0 : set.d.norm();
        lap("trivial_map", t_start);
        return {M, rep};
    }

    // ---- SparseQR (COLAMD) -------------------------------------------------
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> sqr;
    sqr.setPivotThreshold(0.0);
    logging::info(true,
    "SparseQR input: rows=", C_use.rows(),
    " cols=", C_use.cols(),
    " nnz=",  C_use.nonZeros());
    sqr.compute(C_use);
    logging::error(sqr.info() == Eigen::Success, "[ConstraintBuilder] SparseQR failed.");
    lap("SparseQR.compute", t_start);

    Perm P = sqr.colsPermutation();
    SparseMatrix Rsp = sqr.matrixR();
    lap("extract_R", t_start);

    // ---- Rank determination -------------------------------------------------
    const int rmax = std::min<int>(Rsp.rows(), Rsp.cols());
    Precision max_abs_diag = 0;
    for (int k = 0; k < rmax; ++k) {
        Precision akk = 0;
        for (SparseMatrix::InnerIterator it(Rsp, k); it; ++it) {
            if (it.row() == k) { akk = std::abs(it.value()); break; }
            if (it.row() >  k)  break;
        }
        if (akk > max_abs_diag) max_abs_diag = akk;
    }
    rep.R11_max_diag = max_abs_diag;
    const Precision tol = std::max<Precision>(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);
    int r = 0;
    for (; r < rmax; ++r) {
        Precision akk = 0;
        for (SparseMatrix::InnerIterator it(Rsp, r); it; ++it) {
            if (it.row() == r) { akk = std::abs(it.value()); break; }
            if (it.row() >  r)  break;
        }
        if (akk <= tol) break;
    }
    rep.rank = r;
    const int nm_use = n_use - r;
    lap("rank_determination", t_start);

    // ---- Build Xd -----------------------------------------------------------
    Mat Xd = Mat::Zero(r, nm_use);
    if (r > 0 && nm_use > 0) {
        Mat R11 = Mat::Zero(r, r);
        Mat R12 = Mat::Zero(r, nm_use);
        const int Rcols = Rsp.cols();
        for (int j = 0; j < Rcols; ++j) {
            for (SparseMatrix::InnerIterator it(Rsp, j); it; ++it) {
                const int i = it.row();
                if (i >= r) continue;
                const Precision v = it.value();
                if (j < r)             R11(i, j)     = v;
                else if (j < r+nm_use) R12(i, j-r)   = v;
            }
        }
        Xd.noalias() = -R11.template triangularView<Eigen::Upper>().solve(R12);

        if (opt.threshold_X) {
            for (int j = 0; j < nm_use; ++j) {
                const Precision cn = Xd.col(j).norm();
                if (cn == Precision(0)) continue;
                const Precision thr = cn * opt.X_drop_tol_rel;
                for (int i = 0; i < r; ++i)
                    if (std::abs(Xd(i,j)) < thr) Xd(i,j) = Precision(0);
            }
        }
    }
    lap("build_Xd", t_start);

    // ---- Partition columns --------------------------------------------------
    const auto& p = P.indices();
    std::vector<int> slaves_loc(r), masters_loc(nm_use);
    for (int i=0;i<r;++i)      slaves_loc[i]  = p(i);
    for (int j=0;j<nm_use;++j) masters_loc[j] = p(r+j);
    lap("partition_columns", t_start);

    // ---- Build T ------------------------------------------------------------
    std::vector<char> is_used(set.n, 0);
    for (int c: used) is_used[c] = 1;
    int n_never = 0; for (int g=0; g<set.n; ++g) if(!is_used[g]) ++n_never;
    const int nm_total = nm_use + n_never;

    std::vector<Index> masters_glob; masters_glob.reserve(nm_total);
    std::vector<Index> slaves_glob;  slaves_glob.reserve(r);
    for (int j=0;j<nm_use;++j) masters_glob.push_back(used[masters_loc[j]]);
    for (int g=0;g<set.n;++g) if(!is_used[g]) masters_glob.push_back(g);
    for (int i=0;i<r;++i) slaves_glob.push_back(used[slaves_loc[i]]);

    size_t nnz_T = n_never + nm_use;
    for (int j=0;j<nm_use;++j)
        for (int i=0;i<r;++i)
            if (Xd(i,j)!=0) ++nnz_T;

    SparseMatrix T_full(set.n, nm_total);
    {
        std::vector<Eigen::Triplet<Precision>> tTrips;
        tTrips.reserve(nnz_T);
        for (int j=0;j<nm_use;++j) {
            tTrips.emplace_back(used[masters_loc[j]], j, Precision(1));
            for (int i=0;i<r;++i) {
                Precision v = Xd(i,j);
                if (v!=0) tTrips.emplace_back(used[slaves_loc[i]], j, v);
            }
        }
        int col_off = nm_use;
        for (int g=0;g<set.n;++g)
            if(!is_used[g]) tTrips.emplace_back(g, col_off++, Precision(1));
        T_full.setFromTriplets(tTrips.begin(), tTrips.end());
        T_full.makeCompressed();
    }
    lap("build_T", t_start);

    // ---- Build X ------------------------------------------------------------
    SparseMatrix X(r, nm_total);
    if (r > 0 && nm_use > 0) {
        size_t nnz_X=0;
        for (int j=0;j<nm_use;++j)
            for (int i=0;i<r;++i)
                if (Xd(i,j)!=0) ++nnz_X;
        std::vector<Eigen::Triplet<Precision>> xTrips;
        xTrips.reserve(nnz_X);
        for (int j=0;j<nm_use;++j)
            for (int i=0;i<r;++i) {
                Precision v=Xd(i,j);
                if(v!=0) xTrips.emplace_back(i,j,v);
            }
        X.setFromTriplets(xTrips.begin(), xTrips.end());
        X.makeCompressed();
    }
    lap("build_X", t_start);

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
    rep.d_norm            = (set.d.size()==0)?0:set.d.norm();
    rep.n_redundant_rows  = (rep.m>=rep.rank)?(rep.m-rep.rank):0;
    lap("populate_map", t_start);

    // ---- Inhomogeneous case -------------------------------------------------
    if (!rep.homogeneous) {
        Vec d_dense(set.m); for (int i=0;i<set.m;++i) d_dense[i]=set.d[i];
        Vec u_use = sqr.solve(d_dense);
        DynamicVector u_any=DynamicVector::Zero(set.n);
        for(int k=0;k<n_use;++k) u_any[used[k]]=u_use[k];
        DynamicVector resid=set.C*u_any-set.d;
        rep.residual_norm=resid.norm();
        const Precision feas_tol=opt.feas_tol_rel*(rep.d_norm>0?rep.d_norm:Precision(1));
        rep.feasible=(rep.residual_norm<=feas_tol);
        if(rep.feasible) {
            DynamicVector q_any(M.nm_);
            for(int j=0;j<M.nm_;++j) q_any[j]=u_any[M.masters_[j]];
            DynamicVector Xq=M.X_*q_any;
            for(int i=0;i<M.r_;++i) {
                const Index gi=M.slaves_[i];
                M.u_p_[gi]=u_any[gi]-Xq[i];
            }
        }
    }
    lap("inhomogeneous_case", t_start);

    return {M,rep};
}

} // namespace fem::constraint
