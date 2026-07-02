#pragma once

#include "../solve_device.h"
#include "../../core/types_eig.h"

namespace fem::solver {

DynamicMatrix solve_indirect(SolverDevice device,
                             SparseMatrix& mat,
                             const DynamicMatrix& rhs);

inline DynamicVector solve_indirect(SolverDevice device,
                                    SparseMatrix& mat,
                                    const DynamicVector& rhs) {
    const DynamicMatrix rhs_matrix = rhs;
    const DynamicMatrix solution = solve_indirect(device, mat, rhs_matrix);
    return solution.col(0);
}

inline DynamicVector solve_iter(SolverDevice device,
                                SparseMatrix& mat,
                                const DynamicVector& rhs) {
    return solve_indirect(device, mat, rhs);
}

inline DynamicMatrix solve_iter(SolverDevice device,
                                SparseMatrix& mat,
                                const DynamicMatrix& rhs) {
    return solve_indirect(device, mat, rhs);
}

namespace detail {

DynamicMatrix solve_indirect_cpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs);

DynamicMatrix solve_indirect_gpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs);

} // namespace detail
} // namespace fem::solver
