/******************************************************************************
 * @file particular_solution.cpp
 * @brief Compute min-norm particular solution and project onto affine map.
 ******************************************************************************/

#include "particular_solution.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/SparseQR>

namespace fem { namespace constraint {

ParticularOutput compute_particular_and_project(const ParticularInput& in,
                                                const Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>& qr,
                                                const ConstraintMap& M,
                                                Precision feas_tol_rel,
                                                Precision d_norm)
{
    ParticularOutput out{};
    const int n     = (int)M.n_;
    const int n_use = (int)in.used.size();

    out.u_p = DynamicVector::Zero(n);

    if (in.homogeneous) {
        // Only fixed DOFs contribute
        for (int j = 0; j < n; ++j)
            if (in.is_fixed_col[j]) out.u_p[j] = in.fixed_val[j];
        out.residual_norm = Precision(0);
        out.feasible      = true;
        return out;
    }

    // Build RHS for QR system: rows of C_use align with original row indexing
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> d_dense(in.C_use.rows());
    for (int i = 0; i < in.C_use.rows(); ++i) d_dense[i] = in.d_mod[i];

    // Min-norm solution in reduced space
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> u_use = qr.solve(d_dense);

    // Lift to full u_any (putting values only on 'used' columns)
    DynamicVector u_any = DynamicVector::Zero(n);
    for (int k = 0; k < n_use; ++k) u_any[in.used[k]] = u_use[k];

    // fixed DOFs override
    for (int j = 0; j < n; ++j) if (in.is_fixed_col[j]) u_any[j] = in.fixed_val[j];

    // Residual wrt *reduced* system (no access to full C here)
    // r = C_use * u_use - d_mod
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> r = in.C_use * u_use - d_dense;
    out.residual_norm = r.norm();

    const Precision feas_tol = feas_tol_rel * (d_norm > 0 ? d_norm : Precision(1));
    out.feasible = (out.residual_norm <= feas_tol);

    // Project onto affine map u = u_p + T q and store only slave offsets in u_p
    // (masters are represented by q; fixed already placed into u_p)
    // q_any = u_any[masters]
    DynamicVector q_any(M.nm_);
    for (int j = 0; j < M.nm_; ++j) q_any[j] = u_any[M.masters_[j]];

    // Xq = X * q_any  (size r), then for each slave i: u_p[slave_i] = u_any[slave_i] - Xq[i]
    DynamicVector Xq = M.X_ * q_any;
    for (int i = 0; i < M.r_; ++i) {
        const Index gi = M.slaves_[i];
        out.u_p[gi] = u_any[gi] - Xq[i];
    }

    // Ensure fixed values are kept
    for (int j = 0; j < n; ++j) if (in.is_fixed_col[j]) out.u_p[j] = in.fixed_val[j];

    return out;
}

}} // namespace fem::constraint
