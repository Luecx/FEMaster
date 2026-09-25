#include "get_solver_name.h"

namespace fem::solver {

std::string get_solver_name(SolverDevice device,
                            SolverMethod method,
                            DirectSolverMatrixType matrix_type) {
#ifndef SUPPORT_GPU
    device = CPU;
#endif

    if (method == INDIRECT) {
        return device == GPU
            ? "GPU PCG (CUDA incomplete Cholesky)"
            : "CPU PCG (Eigen IncompleteCholesky)";
    }

    if (device == GPU) {
#ifdef USE_CUDSS
        if (matrix_type == DirectSolverMatrixType::SPD) {
            return "GPU DIRECT cuDSS Cholesky";
        }
        if (matrix_type == DirectSolverMatrixType::Symmetric) {
            return "GPU DIRECT cuDSS LDLT";
        }
        return "GPU DIRECT cuDSS LU";
#else
        if (matrix_type == DirectSolverMatrixType::SPD) {
            return "GPU DIRECT cuSolver Cholesky";
        }
        // The legacy GPU direct path cannot factor symmetric-indefinite or
        // general matrices; solve_direct_gpu falls back to the CPU backend.
        device = CPU;
#endif
    }

#ifdef USE_MKL
    return matrix_type == DirectSolverMatrixType::General
        ? "CPU DIRECT MKL PardisoLU"
        : "CPU DIRECT MKL PardisoLDLT";
#elif defined(USE_ACCELERATE)
    return matrix_type == DirectSolverMatrixType::General
        ? "CPU DIRECT Eigen SparseLU"
        : "CPU DIRECT Apple Accelerate LDLT";
#else
    return matrix_type == DirectSolverMatrixType::General
        ? "CPU DIRECT Eigen SparseLU"
        : "CPU DIRECT Eigen SimplicialLDLT";
#endif
}

} // namespace fem::solver
