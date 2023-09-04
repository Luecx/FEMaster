#pragma once

#include "../cuda/cuda_csr.h"
#include "../cuda/cuda_array.h"
#include "../cuda/cuda_vec.h"
#include "../core/core.h"

#include <string>

namespace solver {

enum SolverDevice{
    GPU,
    CPU
};

enum SolverMethod{
    DIRECT,
    INDIRECT
};

//linalg::DenseNMatrix<float> solve(array::Device device,
//                                  linalg::DenseNMatrix<float>& mat,
//                                  linalg::DenseNMatrix<float>& rhs);
//
DynamicVector solve(SolverDevice device,
                    SparseMatrix& mat,
                    DynamicVector& rhs);

DynamicVector solve_iter(SolverDevice device,
                         SparseMatrix& mat,
                         DynamicVector& rhs);

}    // namespace solver
