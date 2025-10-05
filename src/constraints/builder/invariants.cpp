/**
 * @file invariants.cpp
 * @brief Implements debug-time checks for constraint map invariants.
 *
 * The checks compare matrix-free operations against explicit matrix products
 * to ensure implementation correctness during development builds.
 *
 * @see src/constraints/builder/invariants.h
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "invariants.h"

#include "../../core/logging.h"

#include <algorithm>

namespace fem {
namespace constraint {

/**
 * @copydoc check_invariants
 */
void check_invariants(const ConstraintMap& map, const SparseMatrix& C) {
    const int nm = static_cast<int>(map.nm_);
    DynamicVector basis = DynamicVector::Zero(nm);
    double max_norm = 0.0;
    for (int k = 0; k < nm; ++k) {
        basis.setZero();
        basis[k] = 1;
        DynamicVector cu = C * (map.T_ * basis);
        max_norm = std::max(max_norm, static_cast<double>(cu.norm()));
    }
    logging::error(max_norm <= 1e-12,
                   "[ConstraintBuilder] invariant failed: max_k ||C*T*e_k|| = ",
                   max_norm);

    {
        DynamicVector q = DynamicVector::Random(map.nm_);
        DynamicVector u_explicit = map.T_ * q;
        DynamicVector u_fast = DynamicVector::Zero(map.n_);
        map.apply_T(q, u_fast);
        logging::error((u_explicit - u_fast).norm() <= 1e-12,
                       "[ConstraintBuilder] apply_T mismatch explicit T.");
    }

    {
        DynamicVector y = DynamicVector::Random(map.n_);
        DynamicVector z_explicit = map.T_.transpose() * y;
        DynamicVector z_fast = DynamicVector::Zero(map.nm_);
        map.apply_Tt(y, z_fast);
        logging::error((z_explicit - z_fast).norm() <= 1e-12,
                       "[ConstraintBuilder] apply_Tt mismatch explicit T^T.");
    }
}

} // namespace constraint
} // namespace fem
