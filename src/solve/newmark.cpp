/**
 * @file newmark.cpp
 * @brief Implementation of implicit Newmark-β with single factorization (with logging, no throws).
 */

#include "newmark.h"

#ifndef USE_MKL
  #include <Eigen/SparseCholesky>   // SimplicialLDLT
#else
  #include <Eigen/PardisoSupport>   // PardisoLDLT
  #include "../core/config.h"      // global_config.max_threads
  #include <mkl.h>
#endif

#include "../core/logging.h"
#include "../core/timer.h"

#include <stdexcept>
#include <cmath>        // std::ceil
#include <limits>
#include <iomanip>      // std::setw, std::setprecision

namespace fem::solver {

//------------------------------------------------------------------------------
// Internal helpers (anonymous namespace)
//------------------------------------------------------------------------------

namespace {

    // Detect whether M is (numerically) diagonal (i.e., lumped mass).
    inline bool is_lumped_mass(const SparseMatrix& M, double eps = 0.0)
    {
        for (int k = 0; k < M.outerSize(); ++k) {
            for (SparseMatrix::InnerIterator it(M, k); it; ++it) {
                const int i = it.row();
                const int j = it.col();
                if (i != j) {
                    if (eps <= 0.0) {
                        if (it.value() != Precision(0)) return false;
                    } else {
                        if (std::abs(it.value()) > eps) return false;
                    }
                }
            }
        }
        return true;
    }

    struct EffectiveSolver {
        using SolverScalar = Precision;

        const SparseMatrix& M;
        const SparseMatrix& C;
        const SparseMatrix& K;

        // Effective matrix A = K + a0*M + a1*C
        SparseMatrix A;

    #ifndef USE_MKL
        Eigen::SimplicialLDLT<SparseMatrix> fact;
    #else
        Eigen::PardisoLDLT<SparseMatrix>     fact;
    #endif

        bool built = false;
        bool ready = false;

        void build(SolverScalar a0, SolverScalar a1) {
            A = K;
            if (a0 != SolverScalar(0)) A = A + a0 * M;
            if (a1 != SolverScalar(0)) A = A + a1 * C;
            built = true;
            ready = false;
        }

        void factorize() {
        #ifdef USE_MKL
            if (global_config.max_threads > 0) {
                mkl_set_num_threads(global_config.max_threads);
            }
        #endif
            fact.compute(A);
            logging::error(fact.info() == Eigen::Success, "Newmark: factorization of effective matrix A failed.");
            ready = true;
        }

        void solve(const DynamicVector& rhs, DynamicVector& x) {
            logging::error(built && ready, "Newmark: EffectiveSolver not built/factorized.");
            x = fact.solve(rhs);
            logging::error(fact.info() == Eigen::Success, "Newmark: linear solve with A failed.");
        }
    };

    // Compute initial acceleration a0: a0 = M^{-1} (f0 - C v0 - K u0)
    static DynamicVector compute_initial_accel(const SparseMatrix& M,
                                               const SparseMatrix& C,
                                               const SparseMatrix& K,
                                               const DynamicVector& u0,
                                               const DynamicVector& v0,
                                               const DynamicVector& f0)
    {
        const auto n = static_cast<Eigen::Index>(u0.size());
        DynamicVector rhs = f0 - C * v0 - K * u0;

        if (is_lumped_mass(M)) {
            logging::info(true, "Initial acceleration: using lumped (diagonal) mass fast path.");
            DynamicVector a0(n);
            a0.setZero();
            const auto Md = M.diagonal();
            for (Eigen::Index i = 0; i < n; ++i) {
                const auto mii = Md[i];
                logging::error(mii != Precision(0), "Newmark: zero diagonal entry in lumped mass.");
                a0[i] = rhs[i] / mii;
            }
            return a0;
        }

        logging::info(true, "Initial acceleration: factorizing consistent mass (one-time).");
        Timer tM{};
        tM.start();
    #ifndef USE_MKL
        Eigen::SimplicialLDLT<SparseMatrix> Mldl;
        Mldl.compute(M);
        logging::error(Mldl.info() == Eigen::Success, "Newmark: factorization of M failed for initial acceleration.");
        auto a0 = Mldl.solve(rhs);
        logging::error(Mldl.info() == Eigen::Success, "Newmark: solve with M failed for initial acceleration.");
    #else
        Eigen::PardisoLDLT<SparseMatrix> Mldl;
        Mldl.compute(M);
        logging::error(Mldl.info() == Eigen::Success, "Newmark: Pardiso factorization of M failed for initial acceleration.");
        auto a0 = Mldl.solve(rhs);
        logging::error(Mldl.info() == Eigen::Success, "Newmark: Pardiso solve with M failed for initial acceleration.");
    #endif
        tM.stop();
        logging::info(true, "Initial acceleration: done in ", tM.elapsed(), " ms.");
        return a0;
    }

} // namespace (anonymous)

//------------------------------------------------------------------------------
// Public API
//------------------------------------------------------------------------------

NewmarkResult
newmark_linear(const SparseMatrix& M,
               const SparseMatrix& C,
               const SparseMatrix& K,
               const NewmarkIC& ic,
               const NewmarkOpts& opts,
               ForceFn f_of_t)
{
    // Basic checks via logging::error (no throws)
    logging::error(opts.dt > 0.0, "Newmark: dt must be positive.");
    logging::error(opts.t_end >= 0.0, "Newmark: t_end must be non-negative.");
    logging::error(ic.u0.size() != 0 && ic.v0.size() != 0, "Newmark: u0 and v0 must be provided.");
    logging::error(ic.u0.size() == ic.v0.size(), "Newmark: u0 and v0 size mismatch.");

    // Device/method fallbacks (GPU/INDIRECT => CPU/DIRECT here)
    (void)opts.device;
    (void)opts.method;

    const double beta  = opts.beta;
    const double gamma = opts.gamma;
    const double dt    = opts.dt;

    // Banner + matrix stats
    const auto N   = static_cast<long long>(M.rows());
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "IMPLICIT NEWMARK-β (fixed Δt) — linear transient");
    logging::info(true, "-----------------------------------------------------------------------------------------------");
    logging::info(true, "n dof  : ", N);
    logging::info(true, "nnz(M) : ", static_cast<long long>(M.nonZeros()));
    logging::info(true, "nnz(C) : ", static_cast<long long>(C.nonZeros()));
    logging::info(true, "nnz(K) : ", static_cast<long long>(K.nonZeros()));
    logging::info(true, "dt     : ", std::scientific, std::setprecision(6), dt);
    logging::info(true, "t_end  : ", std::scientific, std::setprecision(6), opts.t_end);
    logging::info(true, "β, γ   : ", std::fixed, std::setprecision(4), beta, ", ", gamma);
    logging::down();

    // Newmark constants (fixed step)
    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);
    const double a2 = 1.0 / (beta * dt);
    const double a3 = 1.0 / (2.0 * beta) - 1.0;
    const double a4 = gamma / beta - 1.0;
    const double a5 = dt * (gamma / (2.0 * beta) - 1.0);

    // Build and factorize effective matrix once
    EffectiveSolver eff{M, C, K};

    logging::info(true, "");
    logging::info(true, "Building effective matrix A = K + a0*M + a1*C ...");
    Timer tBuild{};
    tBuild.start();
    eff.build(a0, a1);
    tBuild.stop();
    logging::info(true, "A built: nnz(A) = ", static_cast<long long>(eff.A.nonZeros()),
                       ", time = ", tBuild.elapsed(), " ms.");

    logging::info(true, "Factorizing A (single factorization) ...");
    Timer tFact{};
    tFact.start();
    eff.factorize();
    tFact.stop();
    logging::info(true, "Factorization done in ", tFact.elapsed(), " ms.");

    // Initialize state
    DynamicVector u = ic.u0;
    DynamicVector v = ic.v0;
    DynamicVector a;

    if (ic.a0.size() == ic.u0.size()) {
        logging::info(true, "Using user-supplied initial acceleration.");
        a = ic.a0;
    } else {
        logging::info(true, "Computing initial acceleration from equilibrium at t = 0 ...");
        const DynamicVector f0 = f_of_t(0.0);
        a = compute_initial_accel(M, C, K, u, v, f0);
    }

    // Time stepping
    const int n_steps = static_cast<int>(std::ceil(opts.t_end / dt));
    const int print_stride = std::max(1, n_steps / 20); // ~20 prints total

    NewmarkResult out;
    out.t.reserve(n_steps + 1);
    out.u.reserve(n_steps + 1);
    out.v.reserve(n_steps + 1);
    out.a.reserve(n_steps + 1);

    // Store initial state
    out.t.push_back(0.0);
    out.u.push_back(u);
    out.v.push_back(v);
    out.a.push_back(a);

    logging::info(true, "");
    logging::info(true, "Starting Newmark time marching (", n_steps, " steps) ...");

    Timer tAll{};
    tAll.start();

    // Work vectors
    DynamicVector rhs(u.size());
    DynamicVector un(u.size()), an(u.size()), vn(u.size());

    Timer tLoop{};
    tLoop.start();

    for (int k = 0; k < n_steps; ++k) {
        const double tnp1 = (k + 1) * dt;

        // r_{n+1} = f_{n+1} + M(...) + C(...)
        rhs = f_of_t(tnp1);
        rhs += M * (a0 * u + a2 * v + a3 * a);
        rhs += C * (a1 * u + a4 * v + a5 * a);

        // Solve A u_{n+1} = rhs  (reuse factorization)
        eff.solve(rhs, un);

        // Recover acceleration and velocity (local ops)
        an = a0 * (un - u) - a2 * v - a3 * a;
        vn = v + dt * ((1.0 - gamma) * a + gamma * an);

        // Commit
        u.swap(un);
        v.swap(vn);
        a.swap(an);

        // Store
        out.t.push_back(tnp1);
        out.u.push_back(u);
        out.v.push_back(v);
        out.a.push_back(a);

        if (((k + 1) % print_stride) == 0 || (k + 1) == n_steps) {
            logging::info(true, "Newmark step ",
                          std::setw(8), (k + 1), "/", n_steps,
                          "  t=", std::scientific, std::setprecision(6), tnp1, " s");
            logging::down();
        }
    }

    tLoop.stop();
    tAll.stop();

    logging::info(true, "Time marching finished.");
    logging::up();
    logging::info(true, "Total wall time         : ", tAll.elapsed(), " ms");
    logging::info(true, "Loop time (steps only)  : ", tLoop.elapsed(), " ms");
    logging::info(true, "Avg per step            : ",
                         (n_steps > 0 ? tLoop.elapsed() / double(n_steps) : 0.0), " ms/step");
    logging::down();

    return out;
}

} // namespace fem::solver
