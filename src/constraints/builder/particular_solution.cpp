/**
 * @file particular_solution.cpp
 * @brief Computes and projects particular solutions for constraint systems.
 *
 * The solver obtains a minimum-norm solution for inhomogeneous constraints and
 * projects it into the affine map `u = u_p + T q`.
 *
 * @see src/constraints/builder/particular_solution.h
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "particular_solution.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

namespace fem {
namespace constraint {

/**
 * @copydoc compute_particular_and_project
 */
ParticularOutput compute_particular_and_project(const ParticularInput& input,
                                                const Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>& qr,
                                                const ConstraintMap& map,
                                                Precision feas_tol_rel,
                                                Precision d_norm) {
    ParticularOutput output;
    const int n = static_cast<int>(map.n_);
    const int n_use = static_cast<int>(input.used.size());

    output.u_p = DynamicVector::Zero(n);

    if (input.homogeneous) {
        for (int j = 0; j < n; ++j) {
            if (input.is_fixed_col[j]) {
                output.u_p[j] = input.fixed_val[j];
            }
        }
        output.residual_norm = Precision(0);
        output.feasible = true;
        return output;
    }

    Eigen::VectorXd d_dense(input.C_use.rows());
    for (int i = 0; i < input.C_use.rows(); ++i) {
        d_dense[i] = input.d_mod[i];
    }

    Eigen::VectorXd u_use = qr.solve(d_dense);

    DynamicVector u_any = DynamicVector::Zero(n);
    for (int k = 0; k < n_use; ++k) {
        u_any[input.used[k]] = u_use[k];
    }

    for (int j = 0; j < n; ++j) {
        if (input.is_fixed_col[j]) {
            u_any[j] = input.fixed_val[j];
        }
    }

    Eigen::VectorXd residual = input.C_use * u_use - d_dense;
    output.residual_norm = residual.norm();

    const Precision feas_tol = feas_tol_rel * (d_norm > 0 ? d_norm : Precision(1));
    output.feasible = (output.residual_norm <= feas_tol);

    DynamicVector q_any(map.nm_);
    for (int j = 0; j < map.nm_; ++j) {
        q_any[j] = u_any[map.masters_[j]];
    }

    DynamicVector Xq = map.X_ * q_any;
    for (int i = 0; i < map.r_; ++i) {
        const Index gi = map.slaves_[i];
        output.u_p[gi] = u_any[gi] - Xq[i];
    }

    for (int j = 0; j < n; ++j) {
        if (input.is_fixed_col[j]) {
            output.u_p[j] = input.fixed_val[j];
        }
    }

    return output;
}

} // namespace constraint
} // namespace fem
