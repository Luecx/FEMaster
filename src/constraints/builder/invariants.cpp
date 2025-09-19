/******************************************************************************
* @file invariants.cpp
 * @brief Optional debug invariants for T and fast-apply functions.
 ******************************************************************************/

#include "invariants.h"

#include "../../core/logging.h"
#include <algorithm>

namespace fem { namespace constraint {

void check_invariants(const ConstraintMap& M, const SparseMatrix& C)
{
    // 1) C*T ≈ 0
    {
        const int nm = (int)M.nm_;
        DynamicVector e = DynamicVector::Zero(nm);
        double max_norm = 0.0;
        for (int k = 0; k < nm; ++k) {
            e.setZero(); e[k] = 1;
            DynamicVector cu = C * (M.T_ * e);
            max_norm = std::max(max_norm, (double)cu.norm());
        }
        logging::error(max_norm <= 1e-12, "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ", max_norm);
    }

    // 2) apply_T vs explicit T
    {
        DynamicVector q = DynamicVector::Random(M.nm_);
        DynamicVector u_explicit = M.T_ * q;
        DynamicVector u_fast     = DynamicVector::Zero(M.n_);
        M.apply_T(q, u_fast);
        logging::error( (u_explicit - u_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_T mismatch explicit T.");
    }

    // 3) apply_Tt vs explicit T^T
    {
        DynamicVector y = DynamicVector::Random(M.n_);
        DynamicVector z_explicit = M.T_.transpose() * y;
        DynamicVector z_fast     = DynamicVector::Zero(M.nm_);
        M.apply_Tt(y, z_fast);
        logging::error( (z_explicit - z_fast).norm() <= 1e-12,
                        "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    }
}

}} // namespace fem::constraint
