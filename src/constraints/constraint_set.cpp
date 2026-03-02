/**
 * @file constraint_set.cpp
 * @brief Implements assembly of the sparse constraint system `C u = d`.
 *
 * The assembly routine converts high-level constraint equations into a sparse
 * matrix and right-hand side suitable for solver consumption.
 *
 * @see src/constraints/constraint_set.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_set.h"

#include "../core/logging.h"

namespace fem {
namespace constraint {

/**
 * @copydoc ConstraintSet::assemble
 */
void ConstraintSet::assemble(const SystemDofIds& dofs, Index n_dofs) {
    dof_map = &dofs;
    n = n_dofs;

    TripletList trips;
    trips.reserve(64 * equations.size());

    kept_row_ids.clear();
    kept_row_ids.reserve(equations.size());

    std::vector<Precision> d_vals;
    d_vals.reserve(equations.size());

    Index row = 0;
    for (Index ei = 0; ei < static_cast<Index>(equations.size()); ++ei) {
        const auto& eq = equations[ei];
        bool has_entry = false;

        for (const auto& term : eq.entries) {
            const int dof_id = dofs(term.node_id, term.dof);
            if (dof_id < 0) {
                continue;
            }
            has_entry = true;
            trips.emplace_back(row, dof_id, term.coeff);
        }

        if (has_entry || opt.zero_row_drop_tol <= Precision(0)) {
            kept_row_ids.push_back(ei);
            d_vals.push_back(eq.rhs);
            ++row;
        }
    }

    m = row;

    C.resize(m, n);
    C.setFromTriplets(trips.begin(), trips.end());
    d = Eigen::Map<DynamicVector>(d_vals.data(), static_cast<Eigen::Index>(d_vals.size()));

    col_scale = DynamicVector::Ones(n);
    row_scale = DynamicVector::Ones(m);

    logging::info(true,
                  "[ConstraintSet] Assembled C: m=", m,
                  " n=", n,
                  " nnz=", static_cast<Index>(C.nonZeros()));
}

} // namespace constraint
} // namespace fem
