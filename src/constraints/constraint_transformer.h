/******************************************************************************
 * @file constraint_transformer.h
 * @brief One-stop façade to build and use the constraint null-space map.
 *
 * -----------------------------------------------------------------------------
 * ## What is this?
 *
 * `ConstraintTransformer` is the convenience wrapper that takes your high-level
 * equations (see `equation.h`), assembles the sparse constraint system
 *
 *      C u = d,
 *
 * and immediately builds the **affine null-space map**
 *
 *      u = u_p + T q,
 *
 * via `ConstraintBuilder`. It then exposes a compact API to:
 *   - assemble reduced operators `A = Tᵀ K T`, `B = Tᵀ K_g T`,
 *   - assemble the reduced RHS `b = Tᵀ (f − K u_p)`,
 *   - apply `T` and `Tᵀ`, recover `u` from `q`, and get reactions `r = K u − f`.
 *
 * This is the single header most solvers will include to **enforce constraints
 * without Lagrange multipliers**, keep SPD in the reduced system, and support
 * static, modal, and buckling analyses.
 *
 * -----------------------------------------------------------------------------
 * ## Minimal usage
 * @code{.cpp}
 * using namespace fem::constraint;
 *
 * // 1) Prepare equations and DOF map
 * Equations eqs;
 * eqs.emplace_back( Equation{ { {nA,0, 1.0} }, 0.0 } );               // u_Ax = 0
 * eqs.emplace_back( Equation{ { {nB,0, 1.0}, {nC,0,-1.0} }, 0.0 } );  // u_Bx - u_Cx = 0
 *
 * // 2) Build transformer
 * ConstraintTransformer::BuildOptions opt;
 * opt.set.scale_columns = true; // recommended
 * ConstraintTransformer CT(eqs, dofmap, n_dofs, opt);
 *
 * // 3) Reduce K u = f → A q = b, solve, and recover u
 * SparseMatrix A = CT.assemble_A(K);
 * DynamicVector b = CT.assemble_b(K, f);
 * Eigen::SimplicialLDLT<SparseMatrix> ldlt; ldlt.compute(A);
 * DynamicVector q = ldlt.solve(b);
 * DynamicVector u = CT.recover_u(q);
 *
 * // Diagnostics
 * if (!CT.feasible()) {
 *     // Inspect CT.report() for residual and suspect rows
 * }
 * @endcode
 *
 * -----------------------------------------------------------------------------
 * ## Notes
 * - Redundant constraints are handled automatically (rank detection).
 * - Infeasible (contradictory) sets are detected; see `report()`.
 * - For very large n, consider matrix-free operators: `opA()` / `opB()` with CG.
 *
 * @see constraint_set.h / constraint_set.cpp   (assemble C and d)
 * @see constraint_builder.h / .cpp             (build T, X, u_p and diagnostics)
 * @see constraint_map.h / .cpp                 (apply/assemble/recover utilities)
 * @date    14.09.2025
 * @author  Finn
 ******************************************************************************/

#pragma once

#include "constraint_set.h"
#include "constraint_builder.h"
#include "constraint_map.h"

namespace fem::constraint {

class ConstraintTransformer {
public:
    /// Combined options passed to the two stages: assemble (ConstraintSet) and build (ConstraintBuilder).
    struct BuildOptions {
        ConstraintSet::Options     set;     ///< Scaling / zero-row policy when assembling C and d.
        ConstraintBuilder::Options builder; ///< Rank/feasibility tolerances; X thresholding; reporting.
    };

    /**
     * @brief Construct from equations and system DOF map and immediately build the map.
     *
     * @param eqs     High-level equations to assemble into `C u = d`.
     * @param dofs    System DOF map (node_id, dof) → global index (>=0) or -1 if inactive.
     * @param n_dofs  Total number of global DOFs (columns in C).
     * @param opt     Build options (assembly scaling, rank/feasibility tolerances).
     */
    ConstraintTransformer(const Equations& eqs,
                          const SystemDofIds& dofs,
                          Index n_dofs,
                          const BuildOptions& opt = {})
    {
        set_.equations = eqs;
        set_.assemble(dofs, n_dofs);
        std::tie(map_, report_) = ConstraintBuilder::build(set_, opt.builder);
    }

    // ---- Access to components & diagnostics ---------------------------------

    /// @return Immutable constraint map (u = u_p + T q).
    const ConstraintMap& map() const { return map_; }

    /// @return Assembled set (C, d) and assembly metadata (sizes, scaling).
    const ConstraintSet& set() const { return set_; }

    /// @return Build diagnostics: rank, feasibility, residuals, suspect rows.
    const ConstraintBuilder::Report& report() const { return report_; }

    // ---- Quick accessors ----------------------------------------------------

    bool  homogeneous() const { return report_.homogeneous; }
    bool  feasible()    const { return report_.feasible; }
    Index rank()        const { return report_.rank; }
    Index n_master()    const { return map_.n_master(); }

    // ---- Reduced assemblies (explicit) --------------------------------------

    /// Assemble `A = Tᵀ K T`.
    SparseMatrix assemble_A(const SparseMatrix& K)  const { return map_.assemble_A(K); }

    /// Assemble `B = Tᵀ K_g T` (e.g., M for modal, K_g for buckling).
    SparseMatrix assemble_B(const SparseMatrix& Kg) const { return map_.assemble_B(Kg); }

    /// Assemble `b = Tᵀ (f − K u_p)`.
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
        return map_.assemble_b(K, f);
    }

    // ---- Matrix-free operators ----------------------------------------------

    /// Matrix-free operator for iterative solves: y = Tᵀ K T x.
    ConstraintMap::OpA opA(const SparseMatrix& K)  const { return ConstraintMap::OpA(map_, K); }

    /// Matrix-free operator for generalized EVP: y = Tᵀ K_g T x.
    ConstraintMap::OpB opB(const SparseMatrix& Kg) const { return ConstraintMap::OpB(map_, Kg); }

    // ---- Applications & utilities -------------------------------------------

    /// Apply `T`: u = u_p + T q.
    void apply_T(const DynamicVector& q, DynamicVector& u) const { map_.apply_T(q, u); }

    /// Apply `Tᵀ`: z = Tᵀ y.
    void apply_Tt(const DynamicVector& y, DynamicVector& z) const { map_.apply_Tt(y, z); }

    /// Recover full vector: u = u_p + T q.
    DynamicVector recover_u(const DynamicVector& q) const { return map_.recover_u(q); }

    /// Full residual / reactions: r = K u − f (lives in range(Cᵀ) when solved correctly).
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const {
        return map_.reactions(K, f, q);
    }

private:
    ConstraintSet             set_;
    ConstraintMap             map_;
    ConstraintBuilder::Report report_;
};

} // namespace fem::constraint
