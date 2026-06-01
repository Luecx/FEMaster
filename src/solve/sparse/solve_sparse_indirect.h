#pragma once

#include "../solve_device.h"
#include "../../core/types_eig.h"

namespace fem::solver {

DynamicVector solve_indirect(SolverDevice device,
                             SparseMatrix& mat,
                             DynamicVector& rhs);

inline DynamicVector solve_iter(SolverDevice device,
                                SparseMatrix& mat,
                                DynamicVector& rhs) {
    return solve_indirect(device, mat, rhs);
}

namespace detail {

DynamicVector solve_indirect_cpu(SparseMatrix& mat,
                                 DynamicVector& rhs);

DynamicVector solve_indirect_gpu(SparseMatrix& mat,
                                 DynamicVector& rhs);

} // namespace detail
} // namespace fem::solver
