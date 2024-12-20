#pragma once

#include "../cuda/cuda_csr.h"
#include "../cuda/cuda_array.h"
#include "device.h"  // Import SolverDevice
#include "method.h"  // Import SolverMethod
#include "../core/types_eig.h"
#include "../core/config.h"


namespace fem::solver {

/******************************************************************************
 * @brief Solves a linear system using a direct solver on the specified device.
 *
 * @param device The computational device to use (CPU or GPU).
 * @param mat The sparse matrix representing the system of equations.
 * @param rhs The right-hand side vector.
 * @return DynamicVector The solution vector.
 ******************************************************************************/
DynamicVector solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           DynamicVector& rhs);

/******************************************************************************
 * @brief Solves a linear system using an iterative solver on the specified device.
 *
 * @param device The computational device to use (CPU or GPU).
 * @param mat The sparse matrix representing the system of equations.
 * @param rhs The right-hand side vector.
 * @return DynamicVector The solution vector.
 ******************************************************************************/
DynamicVector solve_iter(SolverDevice device,
                         SparseMatrix& mat,
                         DynamicVector& rhs);

/******************************************************************************
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
 *          INDIRECT, it uses `solve_iter()`.
 ******************************************************************************/
inline DynamicVector solve(SolverDevice device,
                    SolverMethod method,
                    SparseMatrix& mat,
                    DynamicVector& rhs) {
    if (method == DIRECT) {
        return solve_direct(device, mat, rhs);
    } else {
        return solve_iter(device, mat, rhs);
    }
}

}  // namespace solver
