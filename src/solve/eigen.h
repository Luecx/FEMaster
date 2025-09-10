#pragma once

#include "../cuda/cuda_csr.h"
#include "../cuda/cuda_array.h"
#include "../core/core.h"
#include "device.h"

#include <utility>
#include <vector>
#include <optional>

namespace fem::solver {

// ---------- Return type ----------
struct EigenValueVectorPair {
    Precision     value;   // λ
    DynamicVector vector;  // x  (size = n)
};

// ---------- Controls ----------
enum class EigenMode {
    Regular,      // no shift: standard solver
    ShiftInvert,  // (A - σB)^(-1)B (or (A-σI)^(-1) for simple)
    Buckling,     // (A - σB)^(-1)A   (requires A ≻ 0)
    Cayley        // optional specialized transform
};

struct EigenOpts {
    EigenMode mode = EigenMode::Regular;
    double    sigma = 0.0;         // used for ShiftInvert/Buckling/Cayley
    int       ncv = -1;            // Krylov subspace size; -1 => auto
    int       maxit = 3000;
    double    tol = 1e-8;
    enum class Sort { LargestAlge, LargestMagn, SmallestAlge, SmallestMagn } sort = Sort::SmallestAlge;
};

// ---------- Simple problem:  A x = λ x ----------
std::vector<EigenValueVectorPair>
    eigs(SolverDevice device,
         const SparseMatrix& A,
         int k,
         const EigenOpts& opts = {});

// ---------- Generalized problem:  A x = λ B x ----------
std::vector<EigenValueVectorPair>
    eigs(SolverDevice device,
         const SparseMatrix& A,
         const SparseMatrix& B,
         int k,
         const EigenOpts& opts);

} // namespace fem::solver
