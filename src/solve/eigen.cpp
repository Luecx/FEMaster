// src/solve/eigen.cpp
#include "eigen.h"
#include "device.h"
#include "../core/logging.h"
#include "../core/timer.h"

#include <algorithm>
#include <random>
#include <cmath>

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>   // Cholesky backends

// --- Spectra (partial, sparse eigenvalues) ---
#include <Spectra/SymEigsSolver.h>          // simple, regular
#include <Spectra/SymEigsShiftSolver.h>     // simple, shift-invert

#include <Spectra/SymGEigsSolver.h>         // generalized, regular (Cholesky/RegularInverse)
#include <Spectra/SymGEigsShiftSolver.h>    // generalized, shift-invert/buckling/cayley

#include <Spectra/MatOp/SparseSymMatProd.h>     // (A v)
#include <Spectra/MatOp/SparseSymShiftSolve.h>  // (A - σ I)^{-1} v  (simple)
#include <Spectra/MatOp/SymShiftInvert.h>       // (A - σ B)^{-1} * {B|A}  (generalized)

#include <Spectra/MatOp/SparseCholesky.h>       // B-Cholesky (GEigsMode::Cholesky)
// #include <Spectra/MatOp/SparseRegularInverse.h>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#include <mkl.h>
#endif

namespace {

// -----------------------------------------------------------------------------
// Utilities (private to this translation unit)
// -----------------------------------------------------------------------------

/******************************************************************************
 * @brief Map our Sort enum to Spectra::SortRule.
 ******************************************************************************/
inline Spectra::SortRule to_rule(fem::solver::EigenOpts::Sort s) {
    using SR   = Spectra::SortRule;
    using Sort = fem::solver::EigenOpts::Sort;
    switch (s) {
        case Sort::LargestAlge:  return SR::LargestAlge;
        case Sort::LargestMagn:  return SR::LargestMagn;
        case Sort::SmallestAlge: return SR::SmallestAlge;
        default:                 return SR::SmallestMagn;
    }
}

/******************************************************************************
 * @brief Choose a Spectra subspace size (ncv) satisfying k < ncv <= n.
 *
 * Heuristic: ncv = min(n, max(k+2, 2k+20)).
 ******************************************************************************/
inline int choose_ncv(int n, int k, int user) {
    if (user > 0) {
        int ncv = std::min(n, user);
        if (ncv <= k) ncv = std::min(n, std::max(k + 2, k + 1));
        return ncv;
    }
    const int lo  = k + 2;
    const int hi  = 2 * k + 20;
    int ncv       = std::min(n, std::max(lo, hi));
    if (ncv <= k) ncv = std::min(n, k + 2);
    return ncv;
}

/******************************************************************************
 * @brief Ensure an Eigen sparse matrix is compressed (CSR-like) for Spectra.
 ******************************************************************************/
inline void ensure_compressed(const fem::SparseMatrix& M_const) {
    auto& M = const_cast<fem::SparseMatrix&>(M_const);
    if (!M.isCompressed()) M.makeCompressed();
}

/******************************************************************************
 * @brief Relative symmetry error ||A - Aᵀ||_F / ||A||_F without densification.
 ******************************************************************************/
inline double symmetry_error_rel(const fem::SparseMatrix& A)
{
    using Sparse = fem::SparseMatrix;
    double num = 0.0, den = 0.0;

    // ||A - Aᵀ||_F
    {
        Sparse AT = A.transpose();
        Sparse D  = A - AT;
        for (int k = 0; k < D.outerSize(); ++k)
            for (Sparse::InnerIterator it(D, k); it; ++it)
                num += it.value() * it.value();
    }
    // ||A||_F
    for (int k = 0; k < A.outerSize(); ++k)
        for (Sparse::InnerIterator it(A, k); it; ++it)
            den += it.value() * it.value();

    if (den == 0.0) den = 1.0;
    return std::sqrt(num / den);
}

/******************************************************************************
 * @brief Backend-aware SPD test: symmetry + Cholesky acceptance.
 *
 * Uses PardisoLDLT (if USE_MKL) or SimplicialLLT otherwise.
 * No logging; pure predicate.
 *
 * @param A            Square sparse matrix to test.
 * @param rel_sym_tol  Relative symmetry tolerance (default 1e-10).
 * @return true if A is (numerically) symmetric and accepted as SPD.
 ******************************************************************************/
inline bool is_spd(const fem::SparseMatrix& A, double rel_sym_tol = 1e-10)
{
    if (A.rows() != A.cols()) return false;
    const double rel = symmetry_error_rel(A);
    if (!(rel <= rel_sym_tol)) return false;

#ifdef USE_MKL
    {
        Eigen::PardisoLDLT<fem::SparseMatrix> solver;
        solver.compute(A);
        if (solver.info() == Eigen::Success) return true;
        // fall back
    }
#endif
    {
        Eigen::SimplicialLLT<fem::SparseMatrix> llt;
        llt.compute(A);
        if (llt.info() == Eigen::Success) return true;
    }
    return false;
}


} // anonymous

namespace fem::solver {

// -----------------------------
// Simple: A x = λ x
// -----------------------------
std::vector<EigenValueVectorPair>
    eigs(SolverDevice device,
         const SparseMatrix& A,
         int k,
         const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "Matrix must be square for simple eigenproblem");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigen, falling back to CPU");
    device = CPU;
#endif

    const int n = static_cast<int>(A.rows());
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < n, "Requesting k >= N is not supported in sparse partial mode");

#ifdef SUPPORT_GPU
    if (device == GPU) {
        logging::info(true, "GPU path not implemented for eigen; falling back to CPU");
        device = CPU;
    }
#endif

    ensure_compressed(A);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{};
    t.start();

    using OpA = Spectra::SparseSymMatProd<Precision>;
    const int ncv = choose_ncv(n, k, opts.ncv);

    if (opts.mode == EigenMode::Regular) {
        logging::info(true, "Using Spectra::SymEigsSolver (Regular), ncv=", ncv);
        OpA opA(A);
        Spectra::SymEigsSolver<OpA> eig(opA, k, ncv);
        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra partial eigen (Regular) failed");

        DynamicVector evals = eig.eigenvalues();
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::ShiftInvert) {
        logging::info(true, "Using Spectra::SymEigsShiftSolver (ShiftInvert, σ=", opts.sigma, "), ncv=", ncv);
        Spectra::SparseSymShiftSolve<Precision> op(A);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<Precision>>
            eig(op, k, ncv, static_cast<Precision>(opts.sigma));

        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra partial eigen (ShiftInvert) failed");

        DynamicVector evals = eig.eigenvalues();   // already back-transformed
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else {
        logging::error(false, "Simple eigenproblem supports only Regular or ShiftInvert modes");
    }

    t.stop();
    logging::info(true, "Eigen computation (simple) finished on CPU; elapsed: ", t.elapsed(), " ms");
    return out;
}

// -------------------------------------------
// Generalized: A x = λ B x
// -------------------------------------------
std::vector<EigenValueVectorPair>
    eigs(SolverDevice device,
         const SparseMatrix& A,
         const SparseMatrix& B,
         int k,
         const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "A must be square");
    logging::error(B.rows() == B.rows(), "B must be square");
    logging::error(A.rows() == B.rows(), "A and B must have same size");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigen, falling back to CPU");
    device = CPU;
#endif

    const int n = static_cast<int>(A.rows());
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < n, "Requesting k >= N is not supported in sparse partial mode");

#ifdef SUPPORT_GPU
    if (device == GPU) {
        logging::info(true, "GPU path not implemented for generalized eigen; falling back to CPU");
        device = CPU;
    }
#endif

    ensure_compressed(A);
    ensure_compressed(B);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{};
    t.start();

    using OpA = Spectra::SparseSymMatProd<Precision>;
    OpA opA(A);

    const int ncv_user = choose_ncv(n, k, opts.ncv);

    if (opts.mode == EigenMode::Regular) {
        // GEigsMode::Cholesky (requires B ≻ 0)
        logging::info(true, "Using Spectra::SymGEigsSolver (Regular/Cholesky), ncv=", ncv_user, " [B must be SPD]");

        // Explicit SPD check for B
        {
            Eigen::SimplicialLLT<SparseMatrix> llt;
            llt.compute(B);
            logging::error(llt.info() == Eigen::Success,
                           "Cholesky factorization of B failed (B must be SPD)");
        }

        Spectra::SparseCholesky<Precision> opB(B);
        Spectra::SymGEigsSolver<
            OpA,
            Spectra::SparseCholesky<Precision>,
            Spectra::GEigsMode::Cholesky
            > eig(opA, opB, k, ncv_user);

        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (Regular/Cholesky) failed");

        DynamicVector evals = eig.eigenvalues();
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::ShiftInvert) {
        // (A - σB)^{-1} B x = ν x, ν = 1/(λ - σ) [B should be SPD for stability]
        logging::info(true, "Using Spectra::SymGEigsShiftSolver (ShiftInvert, σ=", opts.sigma, "), ncv=", ncv_user,
                      " [B should be SPD]");

        Spectra::SymShiftInvert<Precision>   op(A, B);
        Spectra::SparseSymMatProd<Precision> opB(B);

        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::ShiftInvert
            > eig(op, opB, k, ncv_user, static_cast<Precision>(opts.sigma));

        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (ShiftInvert) failed");

        DynamicVector evals = eig.eigenvalues();   // original λ
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::Buckling) {
        // (A - σB)^{-1} A x = ν x, ν = λ/(λ - σ); Spectra returns physical λ
        // Requirements: A should be SPD; B only symmetric (may be indefinite)
        logging::info(true, "Using Spectra::SymGEigsShiftSolver (Buckling), k=", k);

        // Ensure A is SPD (buckling assumes A ≻ 0)
        logging::error(is_spd(A, 1e-10), "Buckling mode requires A to be SPD");

        // Use a slightly larger subspace for clustered spectra
        const int ncv_buck = std::max(ncv_user, std::min(n, 4 * k + 20));
        logging::info(true, "Buckling: using σ=", opts.sigma, ", ncv=", ncv_buck);

        Spectra::SymShiftInvert<Precision>   op(A, B);   // (A - σB)^{-1} * {...}
        Spectra::SparseSymMatProd<Precision> opA(A);     // IMPORTANT: A-op(K) in buckling

        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Buckling
            > eig(op, opA, k, ncv_buck, static_cast<Precision>(opts.sigma));

        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (Buckling) failed");

        DynamicVector evals = eig.eigenvalues();   // physical λ (buckling factors)
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::Cayley) {
        // Cayley transform via (A - σB)^{-1}{…}
        logging::info(true, "Using Spectra::SymGEigsShiftSolver (Cayley, σ=", opts.sigma, "), ncv=", ncv_user);

        Spectra::SymShiftInvert<Precision>   op(A, B);
        Spectra::SparseSymMatProd<Precision> opB(B);

        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Cayley
            > eig(op, opB, k, ncv_user, static_cast<Precision>(opts.sigma));

        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (Cayley) failed");

        DynamicVector evals = eig.eigenvalues();
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else {
        logging::error(false, "Unsupported EigenMode for generalized eigenproblem");
    }

    t.stop();
    logging::info(true, "Eigen computation (generalized) finished on CPU; elapsed: ", t.elapsed(), " ms");
    return out;
}

} // namespace fem::solver
