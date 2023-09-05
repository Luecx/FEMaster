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

DynamicVector solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           DynamicVector& rhs);

DynamicVector solve_iter(SolverDevice device,
                         SparseMatrix& mat,
                         DynamicVector& rhs);

inline DynamicVector solve(SolverDevice device,
                    SolverMethod method,
                    SparseMatrix& mat,
                    DynamicVector& rhs){
    if (method == DIRECT){
        return solve_direct(device, mat, rhs);
    }else{
        return solve_iter(device, mat, rhs);
    }
}

}    // namespace solver
