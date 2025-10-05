/**
 * @file particular_solution.h
 * @brief Declares routines to compute particular solutions for constraints.
 *
 * The implementation finds a minimum-norm particular solution for inhomogeneous
 * constraints and projects it into the null-space representation.
 *
 * @see src/constraints/builder/particular_solution.cpp
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"
#include "../constraint_map.h"

#include <Eigen/SparseQR>
#include <vector>

namespace fem {
namespace constraint {

/**
 * @struct ParticularInput
 * @brief Input data required to compute the particular solution.
 */
struct ParticularInput {
    bool homogeneous = true;             ///< Whether the system is homogeneous.
    SparseMatrix C_use;                  ///< Reduced constraint matrix.
    DynamicVector d_mod;                 ///< Modified right-hand side aligned to original rows.
    std::vector<int> used;               ///< Mapping from reduced columns to global indices.
    std::vector<char> is_fixed_col;      ///< Flags marking fixed DOFs.
    std::vector<Precision> fixed_val;    ///< Values prescribed on fixed DOFs.
    DynamicVector d_original;            ///< Original right-hand side for residual estimation.
};

/**
 * @struct ParticularOutput
 * @brief Result of the particular-solution computation.
 */
struct ParticularOutput {
    DynamicVector u_p; ///< Particular solution vector of length `n`.
    Precision residual_norm = 0; ///< Norm of the residual `||C u_p - d||`.
    bool feasible = true;        ///< Indicates whether feasibility tolerances are met.
};

/**
 * @brief Computes and projects the particular solution for the constraint system.
 *
 * @param input Particular-solution input data.
 * @param qr Sparse QR factorisation used for solving.
 * @param map Constraint map describing the null space.
 * @param feas_tol_rel Relative feasibility tolerance.
 * @param d_norm Norm of the original right-hand side.
 * @return ParticularOutput Computed particular solution and diagnostics.
 */
ParticularOutput compute_particular_and_project(const ParticularInput& input,
                                                const Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>& qr,
                                                const ConstraintMap& map,
                                                Precision feas_tol_rel,
                                                Precision d_norm);

} // namespace constraint
} // namespace fem
