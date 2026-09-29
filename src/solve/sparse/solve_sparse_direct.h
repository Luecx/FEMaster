#pragma once

#include "../solve_device.h"
#include "../../core/types_eig.h"

namespace fem::solver {

enum class DirectSolverMatrixType {
    SPD,
    Symmetric,
    General
};

struct DirectSolveTimings {
    Time factorization_ms = 0;
    Time backsolve_ms      = 0;
    Time residual_ms       = 0;

    [[nodiscard]] Time total() const {
        return factorization_ms + backsolve_ms + residual_ms;
    }
};

DynamicMatrix solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           const DynamicMatrix& rhs,
                           DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD,
                           DirectSolveTimings* timings = nullptr);

inline DynamicVector solve_direct(SolverDevice device,
                                  SparseMatrix& mat,
                                  const DynamicVector& rhs,
                                  DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD,
                                  DirectSolveTimings* timings = nullptr) {
    const DynamicMatrix rhs_matrix = rhs;
    const DynamicMatrix solution = solve_direct(device, mat, rhs_matrix, matrix_type, timings);
    return solution.col(0);
}

namespace detail {

DynamicMatrix solve_direct_cpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type,
                               DirectSolveTimings* timings = nullptr);

DynamicMatrix solve_direct_gpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type,
                               DirectSolveTimings* timings = nullptr);

} // namespace detail
} // namespace fem::solver
