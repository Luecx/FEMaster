/**
 * @file linear_harmonic.h
 * @brief Declares full and modal-reduced steady-state harmonic response.
 *
 * The load case solves the complex equilibrium equation
 *
 *     (K - omega^2 M + i omega C) u_hat = f_hat,
 *
 * with Rayleigh damping C = alpha M + beta K. Homogeneous support constraints
 * are eliminated before the frequency sweep. A real block representation
 * allows FEMaster's existing real-valued solvers to handle complex responses.
 *
 * The optional *MODALBASIS setting constructs a Ritz space from low-frequency
 * eigenmodes and static residual vectors for the target load patterns.
 * Both approaches recover real and imaginary displacement result fields;
 * corresponding stresses and strains are recovered by the output handler.
 *
 * @see LinearHarmonic
 * @see src/loadcase/linear_harmonic.cpp
 * @author Finn Eggers
 * @date 05.08.2026
 */

#pragma once

#include "loadcase.h"
#include "tools/rayleigh_damping.h"
#include "../constraints/transformer/constraint_transformer.h"
#include "../solve/sparse/solve_sparse.h"

#include <string>
#include <vector>

namespace fem {
namespace loadcase {

/**
 * @brief Executes a frequency sweep for linear structural harmonic response.
 *
 * The default method solves the full constrained complex dynamic stiffness
 * through a real sparse block system. When modal_basis is enabled, low-frequency
 * generalized eigenvectors and mass-orthogonal static load responses define
 * a reduced Ritz space for the same frequency-dependent equations.
 *
 * The step-owned target_conditions supply supports and frequency-dependent
 * loads. Structural matrices and the null-space constraint transformation
 * are reused throughout the sweep. Real and imaginary fields are written
 * separately.
 */
struct LinearHarmonic : public LoadCase {
    // Excitation frequencies in Hz; support and load state is inherited
    // from LoadCase through target_conditions.
    // The solver evaluates one result frame per value
    std::vector<Precision> frequencies;

    // Solver backend for the full or reduced harmonic system
    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;

    // Homogeneous supports are reduced by the null-space transformation
    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    // Rayleigh damping coefficients and optional modal residual reduction
    tools::RayleighDamping damping;
    bool modal_basis = false;


    void set_damping(const tools::RayleighDamping& value) { damping = value; }

    // Analysis identity and execution
    std::string type_name() const override { return "LINEARHARMONIC"; }
    void run() override;
};

} // namespace loadcase
} // namespace fem
