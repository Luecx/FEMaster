#pragma once
/**
 * @file cb_particular.h
 * @brief Compute the affine offset u_p for inhomogeneous constraints (C u = d).
 *
 * Context
 * -------
 * After building the homogeneous map u = u_p + T q with slaves/masters determined
 * by QR, we still need the particular solution u_p that satisfies the inhomogeneous
 * right-hand side. We:
 *   1) Solve a min-norm system on the compressed matrix C_use with SparseQR.
 *   2) Lift the solution back to the full vector u_any (including fixed DOFs).
 *   3) If feasible (‖C u_any − d‖ small), convert u_any into u_p by subtracting
 *      the homogeneous contribution on slave entries: u_p[slave] = u_any[slave] − (X q_any),
 *      where q_any is u_any’s master block.
 *
 * Important
 * ---------
 * - We assume that the ordering of C_use rows matches the original C rows for the
 *   kept rows; i.e., d_kept[i] = d_original[i] is valid.
 * - Fixed columns (from single-NNZ equations) are copied directly into u_p.
 * - If d is zero (homogeneous), u_p is zero except for fixed columns.
 */

#include <Eigen/SparseQR>
#include <vector>
#include "../../core/types_eig.h"

namespace fem { namespace constraint {

struct ParticularOut {
    DynamicVector u_p;        ///< size n (full DOF vector of offsets)
    bool          feasible{true};
    Precision     residual_norm{0};
};

/**
 * @brief Compute inhomogeneous offset u_p using the already-factorized SparseQR.
 *
 * @param C_use        The compressed constraint matrix used in the QR factorization (m'×n_use).
 * @param sqr          SparseQR factor of C_use (already computed).
 * @param C_original   Original constraint matrix C (m×n).
 * @param d_original   Original RHS d (length m).
 * @param used         Mapping from local C_use column index → global DOF index.
 * @param is_fixed_col Flags of size n: 1 if DOF was fixed by single-NNZ row, else 0.
 * @param fixed_val    Values for fixed DOFs (same length as is_fixed_col).
 * @param n            Total number of global DOFs.
 * @param feas_tol_rel Relative feasibility tolerance: ‖C u − d‖ ≤ feas_tol_rel·‖d‖.
 * @param d_norm       Precomputed ‖d‖ (0 if d empty).
 * @param X            r×nm_total coupling from masters to slaves (left block = Xd).
 * @param masters      Global master indices (order matches columns of T/X).
 * @param slaves       Global slave indices (order matches rows 0..r-1 of X).
 *
 * @return ParticularOut with u_p, feasibility flag, and residual norm.
 */
ParticularOut compute_particular(const SparseMatrix& C_use,
                                 const Eigen::SparseQR<SparseMatrix, Eigen::AMDOrdering<int>>& sqr,
                                 const SparseMatrix& C_original,
                                 const DynamicVector& d_original,
                                 const std::vector<int>& used,
                                 const std::vector<char>& is_fixed_col,
                                 const std::vector<Precision>& fixed_val,
                                 int n,
                                 Precision feas_tol_rel,
                                 Precision d_norm,
                                 const SparseMatrix& X,
                                 const std::vector<Index>& masters,
                                 const std::vector<Index>& slaves);

}} // namespace fem::constraint
