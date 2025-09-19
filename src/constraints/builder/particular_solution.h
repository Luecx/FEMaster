/******************************************************************************
* @file particular_solution.h
 * @brief Compute min-norm particular solution u_p for inhomogeneous constraints
 *        and project it onto the affine map u = u_p + T q (store only slave offsets).
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"
#include "../constraint_map.h"
#include <vector>
#include <Eigen/SparseQR>

namespace fem { namespace constraint {

struct ParticularInput {
    bool homogeneous = true;
    SparseMatrix        C_use;
    DynamicVector       d_mod;      ///< aligned to original rows
    std::vector<int>    used;       ///< used→global
    std::vector<char>   is_fixed_col;
    std::vector<Precision> fixed_val;
    DynamicVector       d_original; ///< set.d (for residual)
};

struct ParticularOutput {
    DynamicVector u_p;     ///< size n (fixed DOFs + slave offsets)
    Precision     residual_norm = 0;
    bool          feasible = true;
};

ParticularOutput compute_particular_and_project(const ParticularInput& in,
                                                const Eigen::SparseQR<SparseMatrix, Eigen::AMDOrdering<int>>& qr,
                                                const ConstraintMap& M,
                                                Precision feas_tol_rel,
                                                Precision d_norm);

}} // namespace
