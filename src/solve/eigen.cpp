/******************************************************************************
 * @file eigen.cpp
 * @brief Sparse symmetric partial eigenvalue solvers (CPU, optional MKL).
 *
 * - Symmetric-only code paths (no unsymmetric LU flows in user code).
 * - If -DUSE_MKL is defined, shift–invert operators use MKL PARDISO (LDLT).
 * - Otherwise, shift–invert uses Spectra’s Eigen-based operators.
 * - ncv (subspace size) is auto-chosen; there is no user override.
 ******************************************************************************/

#include "eigen.h"
#include "../core/logging.h"
#include "../core/timer.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>   // SimplicialLLT (for SPD checks)

#include <Spectra/SymEigsSolver.h>          // simple, regular
#include <Spectra/SymEigsShiftSolver.h>     // simple, shift-invert

#include <Spectra/SymGEigsSolver.h>         // generalized, regular (Cholesky mode)
#include <Spectra/SymGEigsShiftSolver.h>    // generalized, shift-invert/buckling/cayley

#include <Spectra/MatOp/SparseSymMatProd.h>     // (A v)
#include <Spectra/MatOp/SparseSymShiftSolve.h>  // (A - σ I)^{-1} v (simple, Eigen backend)
#include <Spectra/MatOp/SparseCholesky.h>       // B-Cholesky for GEigsMode::Cholesky
#include <Spectra/MatOp/SymShiftInvert.h>       // (A - σ B)^{-1}{B|A} (Eigen backend)

#ifdef USE_MKL
  #include <Eigen/PardisoSupport>  // PardisoLDLT
  #include <mkl.h>
#endif

// ============================================================================
// Internal utilities (unnamed namespace)
// ============================================================================
namespace {

/** Map high-level sorting rule to Spectra’s SortRule. */
inline Spectra::SortRule to_rule(fem::solver::EigenOpts::Sort s) {
    using SR   = Spectra::SortRule;
    using Sort = fem::solver::EigenOpts::Sort;
    switch (s) {
        case Sort::LargestAlge:  return SR::LargestAlge;
        case Sort::LargestMagn:  return SR::LargestMagn;
        case Sort::SmallestAlge: return SR::SmallestAlge;
        case Sort::SmallestMagn: default: return SR::SmallestMagn;
    }
}

/** Choose subspace size ncv automatically: k < ncv ≤ n. Heuristic: ncv = min(n, max(k+2, 2k+20)). */
inline int choose_ncv(int n, int k) {
    const int lo  = k + 2;
    const int hi  = 2 * k + 20;
    int ncv       = std::min(n, std::max(lo, hi));
    if (ncv <= k) ncv = std::min(n, k + 2);
    return ncv;
}

/** Ensure an Eigen sparse matrix is in compressed storage (CSC). */
inline void ensure_compressed(const fem::SparseMatrix& M_const) {
    auto& M = const_cast<fem::SparseMatrix&>(M_const);
    if (!M.isCompressed()) M.makeCompressed();
}

} // unnamed namespace

// ============================================================================
// MKL-backed shift–solve operators (symmetric-only, LDLT)
// ============================================================================
#ifdef USE_MKL
namespace mkl_ops {

using fem::Precision;
using fem::SparseMatrix;

/**
 * @brief (A - σ I)^{-1} * x using MKL PARDISO (LDLT). Symmetric-only.
 *
 * Spectra expects:
 *   - using Scalar = ...
 *   - rows(), cols()
 *   - perform_op(const Scalar* x, Scalar* y)
 *   - set_shift(const Scalar& sigma)
 */
struct MKLShiftSolveSimple {
    using Scalar = Precision;  // required by Spectra

    const SparseMatrix& A;
    Precision sigma{0};
    int n{0};
    Eigen::PardisoLDLT<SparseMatrix> ldl;

    MKLShiftSolveSimple(const SparseMatrix& Ain, Precision s)
        : A(Ain), sigma(s), n(static_cast<int>(Ain.rows()))
    {
        refactor_();
    }

    void set_shift(const Precision& s) {
        if (s == sigma) return;
        sigma = s;
        refactor_();
    }

    int rows() const { return n; }
    int cols() const { return n; }

    void perform_op(const Precision* x_in, Precision* y_out) const {
        Eigen::Map<const Eigen::Matrix<Precision, Eigen::Dynamic, 1>> x(x_in, n);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> y = ldl.solve(x);
        std::memcpy(y_out, y.data(), sizeof(Precision) * static_cast<size_t>(n));
    }

private:
    void refactor_() {
        SparseMatrix M = A;                 // M ← A - σ I
        if (!M.isCompressed()) M.makeCompressed();
        std::vector<Eigen::Triplet<Precision>> diag;
        diag.reserve(n);
        for (int i = 0; i < n; ++i) diag.emplace_back(i, i, -sigma);
        SparseMatrix D(n, n);
        D.setFromTriplets(diag.begin(), diag.end());
        M += D;
        M.makeCompressed();

        ldl.compute(M);
        if (ldl.info() != Eigen::Success) {
            throw std::runtime_error("MKLShiftSolveSimple: PARDISO LDLT factorization failed");
        }
    }
};

enum class Right { B, A };

/**
 * @brief (A - σ B)^{-1} * (RightMat * x), with RightMat ∈ {B, A}. Symmetric-only.
 *
 * Also provides set_shift() so Spectra can (re)apply the same/updated σ.
 */
template <Right Rm>
struct MKLShiftSolveGeneral {
    using Scalar = Precision;  // required by Spectra

    const SparseMatrix& A;
    const SparseMatrix& B;
    Precision sigma{0};
    int n{0};
    Eigen::PardisoLDLT<SparseMatrix> ldl;

    MKLShiftSolveGeneral(const SparseMatrix& Ain, const SparseMatrix& Bin, Precision s)
        : A(Ain), B(Bin), sigma(s), n(static_cast<int>(Ain.rows()))
    {
        refactor_();
    }

    void set_shift(const Precision& s) {
        if (s == sigma) return;
        sigma = s;
        refactor_();
    }

    int rows() const { return n; }
    int cols() const { return n; }

    void perform_op(const Precision* x_in, Precision* y_out) const {
        Eigen::Map<const Eigen::Matrix<Precision, Eigen::Dynamic, 1>> x(x_in, n);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> rhs =
            (Rm == Right::B) ? (B * x) : (A * x);
        Eigen::Matrix<Precision, Eigen::Dynamic, 1> y = ldl.solve(rhs);
        std::memcpy(y_out, y.data(), sizeof(Precision) * static_cast<size_t>(n));
    }

private:
    void refactor_() {
        SparseMatrix M = A - sigma * B;   // symmetric by assumption
        if (!M.isCompressed()) M.makeCompressed();

        ldl.compute(M);
        if (ldl.info() != Eigen::Success) {
            throw std::runtime_error("MKLShiftSolveGeneral: PARDISO LDLT factorization failed");
        }
    }
};

} // namespace mkl_ops
#endif // USE_MKL

// ============================================================================
// Public API
// ============================================================================
namespace fem::solver {

// --------------------------------------------------------------------------
// Simple: A x = λ x (A symmetric)
// --------------------------------------------------------------------------
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     int k,
     const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "Matrix must be square");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    // GPU not implemented here — fall back to CPU.
    if (device != CPU) {
        logging::info(true, "Eigen: GPU path not implemented; falling back to CPU");
    }

    ensure_compressed(A);

    const int n   = static_cast<int>(A.rows());
    const int ncv = choose_ncv(n, k);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    using OpA = Spectra::SparseSymMatProd<Precision>;

    if (opts.mode == EigenMode::Regular) {
        logging::info(true, "Eigen (simple): Regular, ncv=", ncv);

        OpA opA(A);
        Spectra::SymEigsSolver<OpA> eig(opA, k, ncv);
        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra simple (Regular) failed");

        DynamicVector evals = eig.eigenvalues();
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::ShiftInvert) {
        logging::info(true, "Eigen (simple): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv);

#ifdef USE_MKL
        mkl_ops::MKLShiftSolveSimple op(A, static_cast<Precision>(opts.sigma));
        Spectra::SymEigsShiftSolver<mkl_ops::MKLShiftSolveSimple>
            eig(op, k, ncv, static_cast<Precision>(opts.sigma));
#else
        // Eigen backend (Simplicial LLT/LDLT underneath). Requires A - σI to factor.
        Spectra::SparseSymShiftSolve<Precision> op(A);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<Precision>>
            eig(op, k, ncv, static_cast<Precision>(opts.sigma));
#endif
        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra simple (ShiftInvert) failed");

        DynamicVector evals = eig.eigenvalues();   // back-transformed to original λ
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else {
        logging::error(false, "Simple eigenproblem supports only Regular or ShiftInvert");
    }

    t.stop();
    logging::info(true, "Eigen (simple) finished; elapsed: ", t.elapsed(), " ms");
    return out;
}

// --------------------------------------------------------------------------
// Generalized: A x = λ B x (A,B symmetric; B typically SPD)
// --------------------------------------------------------------------------
std::vector<EigenValueVectorPair>
eigs(SolverDevice device,
     const SparseMatrix& A,
     const SparseMatrix& B,
     int k,
     const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "A must be square");
    logging::error(B.rows() == B.cols(), "B must be square");
    logging::error(A.rows() == B.rows(), "A and B must have the same size");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    if (device != CPU) {
        logging::info(true, "Eigen (generalized): GPU path not implemented; falling back to CPU");
    }

    ensure_compressed(A);
    ensure_compressed(B);

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = choose_ncv(n, k);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    using OpA = Spectra::SparseSymMatProd<Precision>;
    OpA opA(A);

    if (opts.mode == EigenMode::Regular) {
        // GEigsMode::Cholesky expects B ≻ 0.
        logging::info(true, "Eigen (gen): Regular/Cholesky, ncv=", ncv_user, " [B must be SPD]");

        // Clear error if B is not SPD.
        {
            Eigen::SimplicialLLT<SparseMatrix> llt;
            llt.compute(B);
            logging::error(llt.info() == Eigen::Success,
                           "Cholesky of B failed (B must be SPD) for Regular generalized mode");
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
        // (A - σB)^{-1} B x = ν x, ν = 1/(λ - σ) → Spectra returns λ.
        logging::info(true, "Eigen (gen): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv_user);

#ifdef USE_MKL
        mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::B> op(A, B, static_cast<Precision>(opts.sigma));
        Spectra::SparseSymMatProd<Precision> opB(B);
        Spectra::SymGEigsShiftSolver<
            mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::B>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::ShiftInvert
        > eig(op, opB, k, ncv_user, static_cast<Precision>(opts.sigma));
#else
        Spectra::SymShiftInvert<Precision>   op(A, B);
        Spectra::SparseSymMatProd<Precision> opB2(B);
        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::ShiftInvert
        > eig(op, opB2, k, ncv_user, static_cast<Precision>(opts.sigma));
#endif
        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (ShiftInvert) failed");

        DynamicVector evals = eig.eigenvalues();
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::Buckling) {
        // (A - σB)^{-1} A x = ν x, ν = λ/(λ - σ); Spectra returns physical λ.
        const int ncv_buck = std::max(ncv_user, std::min(n, 4 * k + 20));
        logging::info(true, "Eigen (gen): Buckling, sigma=", opts.sigma, ", ncv=", ncv_buck);

#ifdef USE_MKL
        mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::A> op(A, B, static_cast<Precision>(opts.sigma));
        Spectra::SparseSymMatProd<Precision> opA2(A);
        Spectra::SymGEigsShiftSolver<
            mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::A>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Buckling
        > eig(op, opA2, k, ncv_buck, static_cast<Precision>(opts.sigma));
#else
        Spectra::SymShiftInvert<Precision>   op(A, B);
        Spectra::SparseSymMatProd<Precision> opA2(A);
        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Buckling
        > eig(op, opA2, k, ncv_buck, static_cast<Precision>(opts.sigma));
#endif
        eig.init();
        eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);
        logging::error(eig.info() == Spectra::CompInfo::Successful,
                       "Spectra generalized (Buckling) failed");

        DynamicVector evals = eig.eigenvalues();   // physical λ
        DynamicMatrix evecs = eig.eigenvectors();
        for (int i = 0; i < evals.size(); ++i)
            out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });
    }
    else if (opts.mode == EigenMode::Cayley) {
        // Cayley: Spectra uses same Op as ShiftInvert with Right=B.
        logging::info(true, "Eigen (gen): Cayley, sigma=", opts.sigma, ", ncv=", ncv_user);

#ifdef USE_MKL
        mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::B> op(A, B, static_cast<Precision>(opts.sigma));
        Spectra::SparseSymMatProd<Precision> opB(B);
        Spectra::SymGEigsShiftSolver<
            mkl_ops::MKLShiftSolveGeneral<mkl_ops::Right::B>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Cayley
        > eig(op, opB, k, ncv_user, static_cast<Precision>(opts.sigma));
#else
        Spectra::SymShiftInvert<Precision>   op(A, B);
        Spectra::SparseSymMatProd<Precision> opB2(B);
        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<Precision>,
            Spectra::SparseSymMatProd<Precision>,
            Spectra::GEigsMode::Cayley
        > eig(op, opB2, k, ncv_user, static_cast<Precision>(opts.sigma));
#endif
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
    logging::info(true, "Eigen (generalized) finished; elapsed: ", t.elapsed(), " ms");
    return out;
}

} // namespace fem::solver
