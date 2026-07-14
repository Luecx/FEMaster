/**
 * @file constraint_map.cpp
 * @brief Implements the null-space transformation and reduction utilities.
 *
 * The map represents the transformation `u = u_p + T q`, where `u` contains
 * the full system DOFs and `q` contains the independent reduced DOFs.
 *
 * @see src/constraints/transformer/constraint_map.h
 * @see src/constraints/transformer/null_space.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_map.h"

#include "../../core/logging.h"

namespace fem {
namespace constraint {

void ConstraintMap::recover(const DynamicVector& reduced, DynamicVector& full) const {
    logging::warning(
        reduced.size() == static_cast<Eigen::Index>(n_master()),
        "[ConstraintMap::recover] reduced.size()=", reduced.size(),
        " n_master=", n_master()
    );

    // Initialize the full solution with the particular solution
    full = particular;

    // Insert the independent master DOFs
    for (Index master_id = 0; master_id < n_master(); ++master_id) {
        full[masters[static_cast<std::size_t>(master_id)]] += reduced[master_id];
    }

    // Recover the dependent slave DOFs
    for (Index master_id = 0; master_id < static_cast<Index>(X.outerSize()); ++master_id) {
        for (SparseMatrix::InnerIterator entry(X, master_id); entry; ++entry) {
            const Index slave_id = static_cast<Index>(entry.row());

            full[slaves[static_cast<std::size_t>(slave_id)]] +=
                entry.value() * reduced[master_id];
        }
    }
}

void ConstraintMap::project(const DynamicVector& full, DynamicVector& reduced) const {
    logging::warning(
        full.size() == static_cast<Eigen::Index>(full_size),
        "[ConstraintMap::project] full.size()=", full.size(),
        " full_size=", full_size
    );

    reduced.setZero(n_master());

    // Project the independent master DOFs
    for (Index master_id = 0; master_id < n_master(); ++master_id) {
        reduced[master_id] += full[masters[static_cast<std::size_t>(master_id)]];
    }

    // Add the projected slave contributions
    for (Index master_id = 0; master_id < static_cast<Index>(X.outerSize()); ++master_id) {
        Precision contribution = Precision(0);

        for (SparseMatrix::InnerIterator entry(X, master_id); entry; ++entry) {
            const Index slave_id = static_cast<Index>(entry.row());

            contribution += entry.value()
                          * full[slaves[static_cast<std::size_t>(slave_id)]];
        }

        reduced[master_id] += contribution;
    }
}

DynamicVector ConstraintMap::recover(const DynamicVector& reduced) const {
    DynamicVector full(full_size);
    recover(reduced, full);
    return full;
}

DynamicVector ConstraintMap::project(const DynamicVector& full) const {
    DynamicVector reduced(n_master());
    project(full, reduced);
    return reduced;
}

SparseMatrix ConstraintMap::reduce_matrix(const SparseMatrix& matrix) const {
    SparseMatrix transformed = matrix * T;
    SparseMatrix reduced     = T.transpose() * transformed;

    reduced.makeCompressed();
    return reduced;
}

DynamicVector ConstraintMap::reduce_rhs(const SparseMatrix&  matrix,
                                        const DynamicVector& rhs) const {
    const DynamicVector particular_force = matrix * particular;
    const DynamicVector effective_rhs     = rhs - particular_force;

    return project(effective_rhs);
}

} // namespace constraint
} // namespace fem

