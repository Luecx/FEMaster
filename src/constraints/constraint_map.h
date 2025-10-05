/******************************************************************************
 * @file constraint_map.h
 * @brief Declares the null-space map used to eliminate constrained DOFs.
 *
 * A `ConstraintMap` encapsulates the transformation `u = u_p + T q` that
 * projects reduced unknowns back to the full system. It also provides helper
 * functions to assemble reduced operators and recover reactions.
 *
 * @see src/constraints/constraint_map.cpp
 * @see src/constraints/builder/builder.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"

#include <vector>

namespace fem {
namespace constraint {

/******************************************************************************
 * @class ConstraintMap
 * @brief Stores the null-space representation `u = u_p + T q`.
 ******************************************************************************/
class ConstraintMap {
public:
    /******************************************************************************
     * @brief Default constructor leaves the map uninitialized.
     ******************************************************************************/
    ConstraintMap() = default;

    /// Returns the full system size (`n`).
    Index n_full() const { return n_; }

    /// Returns the number of master DOFs (`n_m`).
    Index n_master() const { return nm_; }

    /// Returns the number of slave DOFs (`r`).
    Index n_slave() const { return r_; }

    /// Provides the indices of master DOFs (length `n_m`).
    const std::vector<Index>& master_idx() const { return masters_; }

    /// Provides the indices of slave DOFs (length `r`).
    const std::vector<Index>& slave_idx() const { return slaves_; }

    /// Accesses the particular solution part `u_p`.
    const DynamicVector& u_p() const { return u_p_; }

    /// Accesses the explicit constraint matrix `X` (when available).
    const SparseMatrix& X() const { return X_; }

    /// Accesses the transformation matrix `T`.
    const SparseMatrix& T() const { return T_; }

    /******************************************************************************
     * @brief Applies the transformation `u = T q + u_p`.
     *
     * @param q Reduced coordinates.
     * @param u Output vector of full coordinates.
     ******************************************************************************/
    void apply_T(const DynamicVector& q, DynamicVector& u) const;

    /******************************************************************************
     * @brief Applies the adjoint transformation `z = T^T y`.
     *
     * @param y Full-space vector.
     * @param z Output vector in the reduced space.
     ******************************************************************************/
    void apply_Tt(const DynamicVector& y, DynamicVector& z) const;

    /******************************************************************************
     * @brief Returns `T q + u_p` without overwriting an existing vector.
     *
     * @param q Reduced coordinates.
     * @return DynamicVector Full coordinate vector.
     ******************************************************************************/
    DynamicVector apply_T(const DynamicVector& q) const;

    /******************************************************************************
     * @brief Returns `T^T y` without overwriting an existing vector.
     *
     * @param y Full-space vector.
     * @return DynamicVector Reduced-space result.
     ******************************************************************************/
    DynamicVector apply_Tt(const DynamicVector& y) const;

    /******************************************************************************
     * @brief Assembles the reduced stiffness matrix `A = T^T K T`.
     *
     * @param K Full-space stiffness matrix.
     * @return SparseMatrix Reduced stiffness matrix.
     ******************************************************************************/
    SparseMatrix assemble_A(const SparseMatrix& K) const;

    /******************************************************************************
     * @brief Assembles the reduced geometric matrix `B = T^T K_g T`.
     *
     * @param Kg Full-space geometric matrix.
     * @return SparseMatrix Reduced geometric matrix.
     ******************************************************************************/
    SparseMatrix assemble_B(const SparseMatrix& Kg) const;

    /******************************************************************************
     * @brief Assembles the reduced right-hand side `b = T^T (f - K u_p)`.
     *
     * @param K Full-space stiffness matrix.
     * @param f Full-space load vector.
     * @return DynamicVector Reduced right-hand side vector.
     ******************************************************************************/
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const;

    /******************************************************************************
     * @class OpA
     * @brief Matrix-free operator that evaluates `y = T^T K T x`.
     ******************************************************************************/
    class OpA {
    public:
        /******************************************************************************
         * @brief Constructs a matrix-free operator for a given stiffness matrix.
         *
         * @param map Constraint map providing `T`.
         * @param K Full-space stiffness matrix.
         ******************************************************************************/
        OpA(const ConstraintMap& map, const SparseMatrix& K);

        /******************************************************************************
         * @brief Performs the matrix-free operation `y = T^T K T x`.
         *
         * @param x Reduced-space input vector.
         * @param y Reduced-space output vector.
         ******************************************************************************/
        void perform_op(const DynamicVector& x, DynamicVector& y) const;

    private:
        const ConstraintMap& map_;
        const SparseMatrix& K_;
        mutable DynamicVector u_full_;
        mutable DynamicVector y_full_;
    };

    /******************************************************************************
     * @class OpB
     * @brief Matrix-free operator that evaluates `y = T^T K_g T x`.
     ******************************************************************************/
    class OpB {
    public:
        /******************************************************************************
         * @brief Constructs a matrix-free operator for a geometric matrix.
         *
         * @param map Constraint map providing `T`.
         * @param Kg Full-space geometric matrix.
         ******************************************************************************/
        OpB(const ConstraintMap& map, const SparseMatrix& Kg);

        /******************************************************************************
         * @brief Performs the matrix-free operation `y = T^T K_g T x`.
         *
         * @param x Reduced-space input vector.
         * @param y Reduced-space output vector.
         ******************************************************************************/
        void perform_op(const DynamicVector& x, DynamicVector& y) const;

    private:
        const ConstraintMap& map_;
        const SparseMatrix& Kg_;
        mutable DynamicVector u_full_;
        mutable DynamicVector y_full_;
    };

    /******************************************************************************
     * @brief Recovers the full solution `u = T q + u_p`.
     *
     * @param q Reduced coordinates.
     * @return DynamicVector Full coordinate vector.
     ******************************************************************************/
    DynamicVector recover_u(const DynamicVector& q) const;

    /******************************************************************************
     * @brief Computes reaction forces associated with constrained DOFs.
     *
     * @param K Full-space stiffness matrix.
     * @param f Full-space load vector.
     * @param q Reduced coordinates.
     * @return DynamicVector Reaction forces in the full space.
     ******************************************************************************/
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const;

public:
    Index n_ = 0;    ///< Total number of DOFs in the original system.
    Index r_ = 0;    ///< Number of constrained (slave) DOFs.
    Index nm_ = 0;   ///< Number of unconstrained (master) DOFs.

    std::vector<Index> masters_; ///< Indices of master DOFs.
    std::vector<Index> slaves_;  ///< Indices of slave DOFs.

    SparseMatrix X_;   ///< Explicit constraint matrix (optional).
    DynamicVector u_p_;///< Particular solution vector.
    SparseMatrix T_;   ///< Null-space basis.
};

} // namespace constraint
} // namespace fem
