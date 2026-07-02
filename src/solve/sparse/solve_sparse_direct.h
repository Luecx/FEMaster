#pragma once

#include "../solve_device.h"
#include "../../core/types_eig.h"

namespace fem::solver {

enum class DirectSolverMatrixType {
    SPD,
    General
};

DynamicMatrix solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           const DynamicMatrix& rhs,
                           DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD);

inline DynamicVector solve_direct(SolverDevice device,
                                  SparseMatrix& mat,
                                  const DynamicVector& rhs,
                                  DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD) {
    const DynamicMatrix rhs_matrix = rhs;
    const DynamicMatrix solution = solve_direct(device, mat, rhs_matrix, matrix_type);
    return solution.col(0);
}

namespace detail {

DynamicMatrix solve_direct_cpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type);

DynamicMatrix solve_direct_gpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type);

} // namespace detail
} // namespace fem::solver
