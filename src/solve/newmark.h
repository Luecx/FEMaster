#pragma once
/**
 * @file newmark.h
 * @brief Linear transient dynamics (implicit Newmark-β) with single factorization.
 *
 * Solves:  M u¨ + C u· + K u = f(t)
 *
 * Assumptions:
 *  - Linear, time-invariant system (K, M, C constant and BC-reduced).
 *  - Fixed time step Δt and fixed (β, γ) → single factorization of
 *        A = K + a0*M + a1*C
 *    where a0 = 1/(β Δt²), a1 = γ/(β Δt).
 *
 * Backends:
 *  - If USE_MKL: PardisoLDLT for factorization/solve.
 *  - Else: Eigen::SimplicialLDLT.
 *
 * Device/Method hooks:
 *  - SolverDevice::CPU supported. GPU currently falls back to CPU.
 *  - SolverMethod::DIRECT recommended (single factorization).
 *  - SolverMethod::INDIRECT currently falls back to DIRECT here.
 *
 * Outputs full time history (u, v, a). For large runs, consider thinning.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 21.10.2025
 */

#include "../core/types_eig.h"  // Precision, DynamicVector/Matrix, SparseMatrix
#include "device.h"             // SolverDevice
#include "method.h"             // SolverMethod

#include <functional>
#include <vector>

namespace fem::solver {

/// Time integrator (extend later with HHT if desired).
enum class TimeIntegrator { NewmarkBeta /*, HHTAlpha*/ };

/// Options for Newmark implicit integration with fixed Δt.
struct NewmarkOpts {
    TimeIntegrator scheme = TimeIntegrator::NewmarkBeta;

    // Average-acceleration (β=1/4, γ=1/2): unconditionally stable, no algorithmic damping.
    double beta  = 0.25;
    double gamma = 0.5;

    double dt    = 1e-3;  ///< Fixed time step (seconds)
    double t_end = 1.0;   ///< End time (seconds)

    // Solver stack control (GPU/INDIRECT => CPU/DIRECT fallback in this implementation)
    SolverDevice device = SolverDevice::CPU;
    SolverMethod method = DIRECT;
};

/// Initial conditions. If a0 not provided (size==0), it will be computed from equilibrium.
struct NewmarkIC {
    DynamicVector u0;  ///< initial displacement
    DynamicVector v0;  ///< initial velocity
    DynamicVector a0;  ///< optional initial acceleration (leave empty to compute)
};

/// Force callback: returns already BC-reduced force vector at time t.
using ForceFn = std::function<DynamicVector(double)>;

/// Result containers (dense; consider thinning in post if memory is a concern).
struct NewmarkResult {
    std::vector<double>        t;
    std::vector<DynamicVector> u;
    std::vector<DynamicVector> v;
    std::vector<DynamicVector> a;
};

/**
 * @brief Implicit Newmark-β transient solver (fixed Δt; single factorization).
 *
 * @param M    Mass matrix (SPD; diagonal or full). Detection of lumped vs. consistent is automatic.
 * @param C    Damping matrix (e.g., Rayleigh). Constant over time.
 * @param K    Stiffness matrix (SPD). Constant over time.
 * @param ic   Initial conditions (a0 computed if empty).
 * @param opts Newmark options (β, γ, dt, t_end, device/method).
 * @param f_of_t Force callback f(t).
 * @return Time history (t, u, v, a) sampled at each step including t=0.
 *
 * @note Re-factorizes the effective matrix only once. If you later support variable Δt/parameters,
 *       you must re-build and re-factorize A accordingly.
 */
NewmarkResult
newmark_linear(const SparseMatrix& M,
               const SparseMatrix& C,
               const SparseMatrix& K,
               const NewmarkIC& ic,
               const NewmarkOpts& opts,
               ForceFn f_of_t);

} // namespace fem::solver
