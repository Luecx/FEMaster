/******************************************************************************
 * @file constraint_map.h
 * @brief Immutable null-space map `u = u_p + T q` with fast application helpers.
 *
 * -----------------------------------------------------------------------------
 * ## What is this?
 *
 * This class represents the **affine null-space parameterization** of all
 * displacements `u ∈ R^n` that satisfy a linear constraint system
 *
 *      C u = d .
 *
 * After running the builder, any admissible solution `u` can be written as
 *
 *      u = u_p + T q,
 *
 * where
 *   - `q ∈ R^(n−r)` are the **independent (master) DOFs**,
 *   - the remaining `r` DOFs are **slaves** expressed via `X` (internal to `T`),
 *   - `u_p` is a **particular solution** if `d ≠ 0` (zero for homogeneous sets).
 *
 * This map lets you reduce a constrained system `K u = f` to an **unconstrained**
 * reduced system
 *
 *      A q = b,   with   A = Tᵀ K T,   b = Tᵀ (f − K u_p),
 *
 * solve for `q`, and recover the full `u = u_p + T q` that satisfies `C u = d`
 * (for feasible constraints, up to floating-point tolerance).
 *
 * -----------------------------------------------------------------------------
 * ## ELI5 intuition
 *
 * Think of some DOFs as **masters** (free knobs you can turn) and the rest as
 * **slaves** (tied by strings to the masters). `T` describes how turning the
 * master knobs moves the whole puppet; `u_p` is the “preload/offset” if the
 * strings have a nonzero target `d`.
 *
 * -----------------------------------------------------------------------------
 * ## Minimal usage example
 * @code{.cpp}
 * using namespace fem::constraint;
 *
 * // Build map via ConstraintBuilder (not shown here):
 * ConstraintMap map = from ConstraintBuilder::build(...).first;
 *
 * // Assemble reduced system for SPD K:
 * SparseMatrix A = map.assemble_A(K);               // A = Tᵀ K T
 * DynamicVector b = map.assemble_b(K, f);           // b = Tᵀ (f − K u_p)
 *
 * // Solve A q = b, then recover u:
 * Eigen::SimplicialLDLT<SparseMatrix> ldlt; ldlt.compute(A);
 * DynamicVector q = ldlt.solve(b);
 * DynamicVector u = map.recover_u(q);               // u = u_p + T q
 *
 * // Optional: compute full residual / reactions
 * DynamicVector r = K * u - f;                      // lives in range(Cᵀ)
 * @endcode
 *
 * -----------------------------------------------------------------------------
 * ## Performance notes
 *
 * - `apply_T(q)` and `apply_Tt(y)` are O(nnz(X) + n) and avoid forming dense data.
 * - `assemble_A/B/b` form reduced objects explicitly; for very large `n`, prefer
 *   **matrix-free** matvecs via the `OpA` / `OpB` nested classes and use CG.
 *
 * -----------------------------------------------------------------------------
 * @see constraint_map.cpp for implementations
 * @see constraint_builder.h/.cpp for construction of this map
 * @date    14.09.2025
 * @author  Finn
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include <vector>

namespace fem::constraint {

/**
 * @class ConstraintMap
 * @brief Immutable affine map `u = u_p + T q` and helpers for reduced assembly.
 *
 * This class is constructed by ConstraintBuilder and then treated as an
 * immutable description of the constraint null-space. It offers:
 *  - fast **application** of `T` and `Tᵀ`,
 *  - **assembly** of reduced operators `A = Tᵀ K T`, `B = Tᵀ K_g T`,
 *  - **RHS** reduction `b = Tᵀ (f − K u_p)`,
 *  - **recovery** of full `u` from reduced `q`,
 *  - a simple **operator** wrapper for matrix-free iterative solves.
 */
class ConstraintMap {
public:
    ConstraintMap() = default;

    // --- Sizes ---------------------------------------------------------------

    /// @return Total number of full DOFs `n`.
    Index n_full() const { return n_; }

    /// @return Number of independent (master) DOFs `n−r`.
    Index n_master() const { return nm_; }

    /// @return Number of slave DOFs `r` (numerical rank of C).
    Index n_slave() const { return r_; }

    // --- Index sets (global DOF indices) ------------------------------------

    /// @return Global DOF indices treated as masters (size = n_master()).
    const std::vector<Index>& master_idx() const { return masters_; }

    /// @return Global DOF indices treated as slaves  (size = n_slave()).
    const std::vector<Index>& slave_idx() const { return slaves_; }

    // --- Map components ------------------------------------------------------

    /// @return Particular solution `u_p` (zero if constraints are homogeneous).
    const DynamicVector& u_p() const { return u_p_; }

    /// @return Sparse slave coefficient matrix `X` of size (r × nm) with rows
    ///         aligned to `slave_idx()`:  `u_slaves = X * q + u_p[slaves]`.
    const SparseMatrix& X() const { return X_; }

    /// @return Full sparse map `T` of size (n × nm). Convenience for explicit assembly.
    const SparseMatrix& T() const { return T_; }

    // --- Application primitives ---------------------------------------------

    /**
     * @brief Apply the map in full form.
     * @param q Reduced unknowns (size nm).
     * @param u Output full vector (size n): `u = u_p + T q`.
     */
    void apply_T(const DynamicVector& q, DynamicVector& u) const;

    /**
     * @brief Apply the transpose map.
     * @param y Full vector (size n).
     * @param z Output reduced vector (size nm): `z = Tᵀ y`.
     */
    void apply_Tt(const DynamicVector& y, DynamicVector& z) const;

    /// @overload Returns the result by value.
    DynamicVector apply_T(const DynamicVector& q) const;

    /// @overload Returns the result by value.
    DynamicVector apply_Tt(const DynamicVector& y) const;

    // --- Reduced operators (assembled) --------------------------------------

    /**
     * @brief Assemble reduced stiffness-like operator.
     * @param K Full matrix (typically SPD).
     * @return `A = Tᵀ K T` (size nm × nm).
     */
    SparseMatrix assemble_A(const SparseMatrix& K) const;

    /**
     * @brief Assemble reduced geometric-stiffness-like operator.
     * @param Kg Full geometric matrix.
     * @return `B = Tᵀ K_g T` (size nm × nm).
     */
    SparseMatrix assemble_B(const SparseMatrix& Kg) const;

    /**
     * @brief Assemble reduced right-hand side.
     * @param K Full matrix K.
     * @param f Full load vector f.
     * @return `b = Tᵀ (f − K u_p)` (size nm).
     */
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const;

    // --- Operator-style matvecs (matrix-free) --------------------------------

    /**
     * @class OpA
     * @brief Matrix-free y = Tᵀ (K (T x)) for iterative solvers (e.g., CG).
     *
     * Internally reuses buffers to avoid repeated allocations.
     */
    class OpA {
    public:
        OpA(const ConstraintMap& map, const SparseMatrix& K);
        /// @brief Perform `y = Tᵀ K T x`.
        void perform_op(const DynamicVector& x, DynamicVector& y) const;
    private:
        const ConstraintMap& map_;
        const SparseMatrix&  K_;
        mutable DynamicVector u_full_, y_full_;
    };

    /**
     * @class OpB
     * @brief Matrix-free y = Tᵀ (K_g (T x)) for eigen/buckling reductions.
     */
    class OpB {
    public:
        OpB(const ConstraintMap& map, const SparseMatrix& Kg);
        /// @brief Perform `y = Tᵀ K_g T x`.
        void perform_op(const DynamicVector& x, DynamicVector& y) const;
    private:
        const ConstraintMap& map_;
        const SparseMatrix&  Kg_;
        mutable DynamicVector u_full_, y_full_;
    };

    // --- Recovery & reactions -----------------------------------------------

    /// @brief Recover full vector `u = u_p + T q`.
    DynamicVector recover_u(const DynamicVector& q) const;

    /**
     * @brief Compute full residual / reaction vector.
     * @param K Full matrix K.
     * @param f Full load vector f.
     * @param q Reduced solution q.
     * @return `r = K u − f` (size n). For a correct reduced solve,
     *         `Tᵀ r ≈ 0` and `r ∈ range(Cᵀ)` (reactions).
     */
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const;

public:
    // ------------------------------------------------------------------------
    // Construction (by ConstraintBuilder only).
    // These members are public by design to keep integration simple and avoid
    // friend-only access friction in downstream code. Treat them as read-only.
    // ------------------------------------------------------------------------
    friend struct ConstraintBuilder;

    Index n_{0};             ///< Total DOFs n.
    Index r_{0};             ///< Number of slaves (rank of C).
    Index nm_{0};            ///< Number of masters (n − r).

    std::vector<Index> masters_; ///< Global DOF ids acting as masters (size nm).
    std::vector<Index> slaves_;  ///< Global DOF ids acting as slaves  (size r).

    SparseMatrix  X_;        ///< (r × nm) slave coefficients (rows aligned with `slaves_`).
    DynamicVector u_p_;      ///< Particular solution (size n), zero for homogeneous sets.
    SparseMatrix  T_;        ///< (n × nm) full sparse map; identity on masters, X on slaves.
};

} // namespace fem::constraint
