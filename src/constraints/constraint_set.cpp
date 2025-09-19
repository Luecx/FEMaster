/******************************************************************************
* @file constraint_set.cpp
 * @brief Assemble the sparse constraint system `C u = d` from high-level equations.
 ******************************************************************************/

#include "./constraint_set.h"

#include "../core/logging.h"

namespace fem::constraint {

void ConstraintSet::assemble(const SystemDofIds& dofs, Index n_dofs)
{
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

        for (const auto& t : eq.entries) {
            const int dof_id = dofs(t.node_id, t.dof);
            if (dof_id < 0) continue;
            has_entry = true;
            trips.emplace_back(row, dof_id, t.coeff);
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

    // Mark identity scaling (no scaling applied)
    col_scale = DynamicVector::Ones(n);
    row_scale = DynamicVector::Ones(m);

    logging::info(true, "[ConstraintSet] Assembled C: m=", m, " n=", n, " nnz=", (Index)C.nonZeros());
}

} // namespace fem::constraint
