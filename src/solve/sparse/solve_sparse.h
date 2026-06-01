#pragma once

#include "../solve_device.h"
#include "../solve_method.h"
#include "solve_sparse_direct.h"
#include "solve_sparse_indirect.h"


namespace fem::solver {

/**
 * @brief Solves a linear system using either a direct or iterative method,
 *        depending on the specified method.
 *
 * @param device The computational device to use (CPU or GPU).
 * @param method The solver method to use (DIRECT or INDIRECT).
 * @param mat The sparse matrix representing the system of equations.
 * @param rhs The right-hand side vector.
 * @return DynamicVector The solution vector.
 *
 * @details If the method is DIRECT, it uses `solve_direct()`. If the method is
 *          INDIRECT, it uses `solve_indirect()`.
 */
inline DynamicVector solve(SolverDevice device,
                    SolverMethod method,
                    SparseMatrix& mat,
                    DynamicVector& rhs,
                    DirectSolverMatrixType direct_matrix_type = DirectSolverMatrixType::SPD) {
    if (method == DIRECT) {
        return solve_direct(device, mat, rhs, direct_matrix_type);
    } else {
        return solve_indirect(device, mat, rhs);
    }
}
}  // namespace solver
