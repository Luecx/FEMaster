/**
 * @file linear_harmonic.h
 * @brief Declares the direct linear harmonic response load case.
 *
 * Solves the steady-state response
 *
 *     (K - omega^2 M + i omega C) u_hat = f_hat
 *
 * over a prescribed set of excitation frequencies. The implementation uses a
 * real 2x2 block system so the existing real sparse solver backends can be
 * reused.
 *
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
 * @struct LinearHarmonic
 * @brief Executes a direct linear harmonic response analysis.
 */
struct LinearHarmonic : public LoadCase {
    std::vector<std::string> supps; ///< Support identifiers applied to the model.
    std::vector<std::string> loads; ///< Harmonic load-amplitude identifiers.
    std::vector<Precision> frequencies; ///< Excitation frequencies.

    solver::SolverDevice device = solver::CPU;
    solver::SolverMethod method = solver::DIRECT;

    constraint::ConstraintTransformer::Method constraint_method =
        constraint::ConstraintTransformer::Method::NullSpace;

    tools::RayleighDamping damping; ///< Proportional viscous damping model.

    void set_damping(const tools::RayleighDamping& value) { damping = value; }

    // Analysis identity and execution
    std::string type_name() const override { return "LINEARHARMONIC"; }
    void run() override;
};

} // namespace loadcase
} // namespace fem
