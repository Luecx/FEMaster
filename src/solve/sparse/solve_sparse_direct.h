#pragma once

#include "../solve_device.h"
#include "../../core/types_eig.h"

namespace fem::solver {

enum class DirectSolverMatrixType {
    SPD,
    General
};

DynamicVector solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           DynamicVector& rhs,
                           DirectSolverMatrixType matrix_type = DirectSolverMatrixType::SPD);

namespace detail {

DynamicVector solve_direct_cpu(SparseMatrix& mat,
                               DynamicVector& rhs,
                               DirectSolverMatrixType matrix_type);

DynamicVector solve_direct_gpu(SparseMatrix& mat,
                               DynamicVector& rhs,
                               DirectSolverMatrixType matrix_type);

} // namespace detail
} // namespace fem::solver
