#include "cb_particular.h"

namespace fem { namespace constraint {

ParticularOut compute_particular(const SparseMatrix& C_use,
                                 const Eigen::SparseQR<SparseMatrix, Eigen::AMDOrdering<int>>& sqr,
                                 const SparseMatrix& C_original,
                                 const DynamicVector& d_original,
                                 const std::vector<int>& used,
                                 const std::vector<char>& is_fixed_col,
                                 const std::vector<Precision>& fixed_val,
                                 int n,
                                 Precision feas_tol_rel,
                                 Precision d_norm,
                                 const SparseMatrix& X,
                                 const std::vector<Index>& masters,
                                 const std::vector<Index>& slaves)
{
    ParticularOut out;
    out.u_p = DynamicVector::Zero(n);

    // Fast path: homogeneous or empty RHS → only fixed DOFs contribute.
    if (d_original.size() == 0 || d_original.lpNorm<Eigen::Infinity>() == Precision(0)) {
        for (int j = 0; j < n; ++j)
            if (is_fixed_col[j]) out.u_p[j] = fixed_val[j];
        out.feasible = true;
        out.residual_norm = 0;
        return out;
    }

    // Build RHS for the kept rows (assumes row alignment with C_use).
    DynamicVector d_kept = DynamicVector::Zero(C_use.rows());
    for (int i = 0; i < C_use.rows(); ++i)
        d_kept[i] = d_original[i];

    // Min-norm solve on C_use
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> u_use = sqr.solve(d_kept);

    // Lift back to full vector u_any
    DynamicVector u_any = DynamicVector::Zero(n);
    for (int k = 0; k < static_cast<int>(used.size()); ++k)
        u_any[used[k]] = u_use[k];
    // Insert fixed DOFs
    for (int j = 0; j < n; ++j)
        if (is_fixed_col[j]) u_any[j] = fixed_val[j];

    // Check feasibility on original system
    DynamicVector resid = C_original * u_any - d_original;
    out.residual_norm = resid.norm();
    const Precision feas_tol = feas_tol_rel * (d_norm > 0 ? d_norm : Precision(1));
    out.feasible = (out.residual_norm <= feas_tol);

    if (out.feasible) {
        // Convert u_any to an affine offset u_p: store only slave offsets beyond X * master part
        // q_any = u_any restricted to masters (in masters’ order)
        DynamicVector q_any(static_cast<int>(masters.size()));
        for (int j = 0; j < static_cast<int>(masters.size()); ++j)
            q_any[j] = u_any[masters[j]];

        DynamicVector Xq = X * q_any; // length r (in the same order as `slaves`)
        for (int i = 0; i < static_cast<int>(slaves.size()); ++i)
            out.u_p[slaves[i]] = u_any[slaves[i]] - Xq[i];

        // Fixed DOFs: keep their values in u_p as well
        for (int j = 0; j < n; ++j)
            if (is_fixed_col[j]) out.u_p[j] = fixed_val[j];
    }

    return out;
}

}} // namespace fem::constraint
