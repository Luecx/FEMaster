/**
 * @file constraint_map.h
 * @brief Declares the transformation between reduced and full constraint coordinates.
 *
 * @see src/constraints/transformer/constraint_map.cpp
 * @see src/constraints/transformer/null_space.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"

#include <vector>

namespace fem {
namespace constraint {

/**
 * @brief Stores the null-space transformation `u = u_p + T q`.
 */
struct ConstraintMap {
    // Full system dimensions
    Index full_size{};            ///< Number of unconstrained system DOFs.

    // Reduced DOF partition
    std::vector<Index> masters{}; ///< Full-system indices of the independent master DOFs.
    std::vector<Index> slaves{};  ///< Full-system indices of the dependent slave DOFs.

    // Null-space transformation
    SparseMatrix  T{};            ///< Transformation from reduced to full coordinates.
    SparseMatrix  X{};            ///< Slave-to-master transformation block.
    DynamicVector particular{};   ///< Particular solution `u_p` of the constraint equations.

    // Reduced system dimensions
    Index n_master() const { return static_cast<Index>(masters.size()); }
    Index n_slave()  const { return static_cast<Index>(slaves.size()); }

    // Coordinate transformations
    void          recover(const DynamicVector& reduced, DynamicVector& full) const;
    void          project(const DynamicVector& full, DynamicVector& reduced) const;
    DynamicVector recover(const DynamicVector& reduced) const;
    DynamicVector project(const DynamicVector& full) const;

    // System reduction
    SparseMatrix  reduce_matrix(const SparseMatrix& matrix) const;
    DynamicVector reduce_rhs(const SparseMatrix& matrix, const DynamicVector& rhs) const;
};

} // namespace constraint
} // namespace fem

