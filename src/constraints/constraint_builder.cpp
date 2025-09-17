/******************************************************************************
 * @file constraint_builder.cpp
 * @brief Build a reduced null-space map u = u_p + T q from C u = d.
 *
 * Component-wise SparseQR (COLAMD, single-thread, no caching):
 *   - Compress zero columns of C → C_use (m x n_use)
 *   - Build column components (Union-Find) where columns share a row
 *   - For each component Ck:
 *       * SparseQR: Ck * Pk = Qk * Rk  (pivotThreshold=0)
 *       * rank rk via diag(Rk) relative threshold
 *       * Xk = -R11^{-1} R12   (thin dense)
 *       * optional per-column drop on Xk
 *       * accumulate T/X triplets with correct global indexing
 *   - Append identities for never-used DOFs as masters
 *   - Inhomogeneous u_p via SparseQR::solve on the full C_use if needed
 *
 * Timings are printed per major step and totals for QR/solve over all components.
 *
 * @date    17.09.2025  (component-wise SparseQR + timings)
 * @author  Finn
 ******************************************************************************/

#include "constraint_builder.h"
#include "constraint_map.h"
#include "../core/logging.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
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

// Simple Union-Find for column components
struct UF {
    std::vector<int> p, r;
    explicit UF(int n=0): p(n), r(n,0) { std::iota(p.begin(), p.end(), 0); }
    int find(int a){ while(p[a]!=a){ p[a]=p[p[a]]; a=p[a]; } return a; }
    void unite(int a,int b){ a=find(a); b=find(b); if(a==b) return; if(r[a]<r[b]) std::swap(a,b); p[b]=a; if(r[a]==r[b]) ++r[a]; }
};

static void build_column_components(const SparseMatrix& C_use,
                                    std::vector<std::vector<int>>& comps)
{
    const int n_use = C_use.cols();
    const int m     = C_use.rows();

    UF uf(n_use);

    // Collect columns per row
    std::vector<std::vector<int>> cols_of_row(m);
    for (int j = 0; j < n_use; ++j) {
        for (SparseMatrix::InnerIterator it(C_use, j); it; ++it) {
            cols_of_row[it.row()].push_back(j);
        }
    }
    // Unite columns that share any row
    for (auto& v : cols_of_row) {
        if (v.size() <= 1) continue;
        const int a0 = v[0];
        for (size_t k = 1; k < v.size(); ++k) uf.unite(a0, v[k]);
    }
    // Bucket columns by root
    std::unordered_map<int, std::vector<int>> buckets;
    buckets.reserve(n_use);
    for (int j = 0; j < n_use; ++j) buckets[uf.find(j)].push_back(j);

    comps.clear();
    comps.reserve(buckets.size());
    for (auto& kv : buckets) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end()); // stable output
        comps.emplace_back(std::move(v));
    }
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

    auto t_global_start = Clock::now();
    auto now = [](){ return Clock::now(); };
    auto ms  = [](Clock::time_point a, Clock::time_point b){
        return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    };

    logging::info(true, "[INFO] Begin of building constraint transformer");
    logging::info(true, "  [ConstraintSet] Assembled C: m=", set.m, " n=", set.n, " nnz=", set.C.nonZeros());

    // ---- 1) Compress columns (only constrained DOFs) -----------------------
    auto t0 = now();
    std::vector<int> used; Eigen::VectorXi old2new; SparseMatrix C_use;
    compress_zero_columns(set.C, used, old2new, C_use);
    auto t1 = now();
    logging::info(true, "  compress_zero_columns : ", ms(t0,t1), " ms");
    logging::info(true, "  SparseQR input (C_use): rows=", C_use.rows(),
                  " cols=", C_use.cols(), " nnz=", C_use.nonZeros());

    ConstraintMap M; Report rep;
    rep.m = set.m; rep.n = set.n;
    rep.homogeneous = (set.d.size()==0) || (set.d.lpNorm<Eigen::Infinity>()==Precision(0));

    const int n_use = (int)used.size();
    if (n_use == 0) {
        M.n_ = set.n; M.r_ = 0; M.nm_ = set.n;
        M.masters_.resize(set.n); std::iota(M.masters_.begin(), M.masters_.end(), 0);
        M.slaves_.clear();
        M.T_.resize(set.n, set.n); M.T_.setIdentity();
        M.X_.resize(0, set.n);
        M.u_p_ = DynamicVector::Zero(set.n);
        rep.rank = 0; rep.n_redundant_rows = rep.m; rep.feasible = true;
        rep.d_norm = (set.d.size()==0)?0:set.d.norm();
        return {M, rep};
    }

    // ---- 2) Column components ----------------------------------------------
    auto t2 = now();
    std::vector<std::vector<int>> comps;
    build_column_components(C_use, comps);
    auto t3 = now();
    logging::info(true, "  build_components (union-find) : ", ms(t2,t3), " ms");
    logging::info(true, "  components: ", (int)comps.size());

    // ---- 3) Per-component QR and accumulation ------------------------------
    auto t4 = now();

    // Global collectors
    std::vector<Index> masters_glob; masters_glob.reserve(n_use);
    std::vector<Index> slaves_glob;  slaves_glob.reserve(std::min(n_use, set.m));

    std::vector<Eigen::Triplet<Precision>> T_trips;
    std::vector<Eigen::Triplet<Precision>> X_trips;
    T_trips.reserve((size_t)C_use.nonZeros()*2);
    X_trips.reserve((size_t)C_use.nonZeros());

    int total_r = 0;
    int total_nm_use = 0;

    // For timings
    long long qr_ms_total = 0;
    long long extract_ms_total = 0;
    long long rank_ms_total = 0;
    long long Xd_ms_total = 0;

    // Process each component
    for (const auto& cols : comps) {
        if (cols.empty()) continue;

        const int nc = (int)cols.size();
        // loc2use: local-comp col -> C_use col
        const std::vector<int> loc2use = cols;
        // Collect rows that touch this component
        std::vector<char> row_flag(C_use.rows(), 0);
        int mr = 0;
        for (int jloc = 0; jloc < nc; ++jloc) {
            int juse = loc2use[jloc];
            for (SparseMatrix::InnerIterator it(C_use, juse); it; ++it) {
                if (!row_flag[it.row()]) { row_flag[it.row()] = 1; ++mr; }
            }
        }
        std::vector<int> rows; rows.reserve(mr);
        for (int i = 0; i < C_use.rows(); ++i) if (row_flag[i]) rows.push_back(i);

        // Build Ck
        SparseMatrix Ck(mr, nc);
        {
            std::vector<Eigen::Triplet<Precision>> trips;
            trips.reserve((size_t)nc * 4);
            std::unordered_map<int,int> row2loc; row2loc.reserve(mr*2);
            for (int r = 0; r < mr; ++r) row2loc[rows[r]] = r;

            for (int jloc = 0; jloc < nc; ++jloc) {
                int juse = loc2use[jloc];
                for (SparseMatrix::InnerIterator it(C_use, juse); it; ++it) {
                    const int iloc = row2loc[it.row()];
                    trips.emplace_back(iloc, jloc, it.value());
                }
            }
            Ck.setFromTriplets(trips.begin(), trips.end());
            Ck.makeCompressed();
        }

        // QR on component
        auto t_qr0 = now();
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> sqr;
        sqr.setPivotThreshold(0.0);
        sqr.compute(Ck);
        auto t_qr1 = now();
        qr_ms_total += ms(t_qr0, t_qr1);
        if (sqr.info() != Eigen::Success) {
            logging::error(false, "[ConstraintBuilder] SparseQR failed on component.");
            continue;
        }

        // Extract Rk
        auto t_ex0 = now();
        SparseMatrix Rk = sqr.matrixR();
        auto t_ex1 = now();
        extract_ms_total += ms(t_ex0, t_ex1);

        // Rank via fast diagonal scan
        auto t_rk0 = now();
        const int rmax = std::min<int>(Rk.rows(), Rk.cols());
        Precision max_abs_diag = 0;
        for (int k = 0; k < rmax; ++k) {
            Precision akk = 0;
            for (SparseMatrix::InnerIterator it(Rk, k); it; ++it) {
                const int i = it.row();
                if (i == k) { akk = std::abs(it.value()); break; }
                if (i >  k)  break;
            }
            if (akk > max_abs_diag) max_abs_diag = akk;
        }
        const Precision tol = std::max<Precision>(Precision(1e-300), max_abs_diag * opt.rank_tol_rel);
        int rk = 0;
        for (; rk < rmax; ++rk) {
            Precision akk = 0;
            for (SparseMatrix::InnerIterator it(Rk, rk); it; ++it) {
                const int i = it.row();
                if (i == rk) { akk = std::abs(it.value()); break; }
                if (i >  rk)  break;
            }
            if (akk <= tol) break;
        }
        auto t_rk1 = now();
        rank_ms_total += ms(t_rk0, t_rk1);

        const int nm_k = nc - rk;

        // Build Xd_k
        Mat Xd_k = Mat::Zero(rk, nm_k);
        auto t_xd0 = now();
        if (rk > 0 && nm_k > 0) {
            Mat R11 = Mat::Zero(rk, rk);
            Mat R12 = Mat::Zero(rk, nm_k);

            const int Rcols = Rk.cols();
            for (int j = 0; j < Rcols; ++j) {
                for (SparseMatrix::InnerIterator it(Rk, j); it; ++it) {
                    const int i = it.row();
                    if (i >= rk) continue;
                    const Precision v = it.value();
                    if (j < rk)             R11(i, j)     = v;
                    else if (j < rk+nm_k)   R12(i, j-rk)  = v;
                }
            }
            Xd_k.noalias() = -R11.template triangularView<Eigen::Upper>().solve(R12);

            if (opt.threshold_X) {
                for (int j = 0; j < nm_k; ++j) {
                    const Precision cn = Xd_k.col(j).norm();
                    if (cn == Precision(0)) continue;
                    const Precision thr = cn * opt.X_drop_tol_rel;
                    for (int i = 0; i < rk; ++i)
                        if (std::abs(Xd_k(i,j)) < thr) Xd_k(i,j) = Precision(0);
                }
            }
        }
        auto t_xd1 = now();
        Xd_ms_total += ms(t_xd0, t_xd1);

        // Permutation partition inside component
        const auto& p = sqr.colsPermutation().indices();
        std::vector<int> slaves_loc_comp;  slaves_loc_comp.reserve(rk);
        std::vector<int> masters_loc_comp; masters_loc_comp.reserve(nm_k);
        for (int i = 0; i < rk;    ++i) slaves_loc_comp.push_back(  p(i)       );
        for (int j = 0; j < nm_k;  ++j) masters_loc_comp.push_back( p(rk + j)  );

        // Global indexing
        const int slave_row_offset  = (int)slaves_glob.size();
        const int master_col_offset = (int)masters_glob.size();

        // Append global slave and master indices (DOF ids)
        for (int i = 0; i < rk; ++i)
            slaves_glob.push_back( used[ loc2use[ slaves_loc_comp[i] ] ] );

        for (int j = 0; j < nm_k; ++j)
            masters_glob.push_back( used[ loc2use[ masters_loc_comp[j] ] ] );

        // T trips for this component (only constrained masters)
        for (int j = 0; j < nm_k; ++j) {
            const int jglob = master_col_offset + j;
            const int row_master_loc = masters_loc_comp[j]; // in [0..nc-1]
            const int row_master_use = loc2use[row_master_loc];
            T_trips.emplace_back( used[row_master_use], jglob, Precision(1) );

            for (int i = 0; i < rk; ++i) {
                const Precision v = Xd_k(i, j);
                if (v == Precision(0)) continue;
                const int row_slave_use = loc2use[ slaves_loc_comp[i] ];
                T_trips.emplace_back( used[row_slave_use], jglob, v );
            }
        }

        // X trips (rows aligned with slaves_glob order)
        for (int j = 0; j < nm_k; ++j) {
            const int jglob = master_col_offset + j;
            for (int i = 0; i < rk; ++i) {
                const Precision v = Xd_k(i, j);
                if (v == Precision(0)) continue;
                X_trips.emplace_back( slave_row_offset + i, jglob, v );
            }
        }

        total_r       += rk;
        total_nm_use  += nm_k;
    }

    auto t5 = now();
    logging::info(true, "  per-component SparseQR total    : ", qr_ms_total,     " ms");
    logging::info(true, "  per-component extract_R total   : ", extract_ms_total," ms");
    logging::info(true, "  per-component rank total        : ", rank_ms_total,   " ms");
    logging::info(true, "  per-component build_Xd total    : ", Xd_ms_total,     " ms");
    logging::info(true, "  accumulate_components           : ", ms(t4,t5),       " ms");

    // ---- 4) Append never-used DOFs as masters ------------------------------
    auto t6 = now();
    std::vector<char> is_used(set.n, 0);
    for (int c : used) is_used[c] = 1;
    int n_never = 0; for (int g = 0; g < set.n; ++g) if (!is_used[g]) ++n_never;

    const int nm_total = total_nm_use + n_never;
    // Append identities for never-used masters (and remember their indices)
    int col_off_never = total_nm_use;
    for (int g = 0; g < set.n; ++g) {
        if (!is_used[g]) {
            masters_glob.push_back(g);
            T_trips.emplace_back(g, col_off_never++, Precision(1));
        }
    }
    auto t7 = now();
    logging::info(true, "  append_never_used               : ", ms(t6,t7), " ms");

    // ---- 5) Finalize T and X -----------------------------------------------
    auto t8 = now();
    SparseMatrix T_full(set.n, nm_total);
    T_full.setFromTriplets(T_trips.begin(), T_trips.end());
    T_full.makeCompressed();

    SparseMatrix X(total_r, nm_total);
    X.setFromTriplets(X_trips.begin(), X_trips.end());
    X.makeCompressed();
    auto t9 = now();
    logging::info(true, "  finalize_T_X                    : ", ms(t8,t9), " ms");

    // ---- 6) Populate map & report ------------------------------------------
    M.n_       = set.n;
    M.r_       = total_r;
    M.nm_      = (Index)masters_glob.size();
    M.masters_ = std::move(masters_glob);
    M.slaves_  = std::move(slaves_glob);
    M.T_       = std::move(T_full);
    M.X_       = std::move(X);
    M.u_p_     = DynamicVector::Zero(set.n);

    rep.rank             = total_r;
    rep.feasible         = true;
    rep.d_norm           = (set.d.size()==0)?0:set.d.norm();
    rep.n_redundant_rows = (rep.m >= rep.rank)?(rep.m - rep.rank):0;

#ifndef NDEBUG
    // Cheap invariants (debug)
    {
        const int nm = (int)M.nm_;
        DynamicVector e(nm), cu; double max_norm = 0.0;
        for (int k = 0; k < nm; ++k) {
            e.setZero(); e[k]=1;
            cu = set.C * (M.T_ * e);
            max_norm = std::max(max_norm, (double)cu.norm());
        }
        logging::error(max_norm <= 1e-12, "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ", max_norm);
    }
    {
        DynamicVector q = DynamicVector::Random(M.nm_);
        DynamicVector u_explicit = M.T_*q + M.u_p_;
        DynamicVector u_fast; M.apply_T(q, u_fast = DynamicVector::Zero(M.n_));
        logging::error( (u_explicit - u_fast).norm() <= 1e-12, "[ConstraintBuilder] apply_T mismatch explicit T.");
    }
    {
        DynamicVector y = DynamicVector::Random(M.n_);
        DynamicVector z_explicit = M.T_.transpose()*y;
        DynamicVector z_fast; M.apply_Tt(y, z_fast = DynamicVector::Zero(M.nm_));
        logging::error( (z_explicit - z_fast).norm() <= 1e-12, "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    }
#endif

    // ---- 7) Inhomogeneous (optional, min-norm) ------------------------------
    auto t10 = now();
    if (!rep.homogeneous) {
        // Fallback: solve once on the full C_use (cheapest to implement)
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> sqr_full;
        sqr_full.setPivotThreshold(0.0);
        sqr_full.compute(C_use);
        Vec d_dense(set.m); for (int i=0;i<set.m;++i) d_dense[i]=set.d[i];
        Vec u_use = sqr_full.solve(d_dense);

        DynamicVector u_any = DynamicVector::Zero(set.n);
        for (int k = 0; k < n_use; ++k) u_any[used[k]] = u_use[k];

        DynamicVector resid = set.C * u_any - set.d;
        rep.residual_norm   = resid.norm();
        const Precision feas_tol = opt.feas_tol_rel * (rep.d_norm > 0 ? rep.d_norm : Precision(1));
        rep.feasible = (rep.residual_norm <= feas_tol);

        if (rep.feasible) {
            DynamicVector q_any(M.nm_);
            for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];
            DynamicVector Xq = M.X_ * q_any;
            for (int i = 0; i < M.r_; ++i) {
                const Index gi = M.slaves_[i];
                M.u_p_[gi] = u_any[gi] - Xq[i];
            }
        }
    }
    auto t11 = now();
    logging::info(true, "  inhomogeneous_case              : ", ms(t10,t11), " ms");

    // ---- Summary ------------------------------------------------------------
    auto t_global_end = Clock::now();
    logging::info(true, "  SUMMARY ------------------------------");
    logging::info(true, "    components              : ", (int)comps.size());
    {
        int max_nc = 0; long long sum_nc = 0;
        for (auto& v : comps) { sum_nc += (int)v.size(); max_nc = std::max(max_nc, (int)v.size()); }
        logging::info(true, "    total constrained cols  : ", (int)sum_nc, "  (n_use=", n_use, ")");
        logging::info(true, "    max component cols      : ", max_nc);
        logging::info(true, "    rank(C)                 : ", rep.rank);
    }
    logging::info(true, "    T nnz                   : ", (long long)M.T_.nonZeros());
    logging::info(true, "    X nnz                   : ", (long long)M.X_.nonZeros());
    logging::info(true, "    QR total                : ", qr_ms_total,     " ms");
    logging::info(true, "    extract_R total         : ", extract_ms_total," ms");
    logging::info(true, "    rank total              : ", rank_ms_total,   " ms");
    logging::info(true, "    build_Xd total          : ", Xd_ms_total,     " ms");
    logging::info(true, "    total build             : ", ms(t_global_start, t_global_end), " ms");

    return {M, rep};
}

} // namespace fem::constraint
