/**
 * @file constraint_transformer.h
 * @brief Declares a facade for building and using constraint maps.
 *
 * The `ConstraintTransformer` bundles the steps required to assemble the
 * constraint set, build the null-space map, and expose convenience helpers for
 * reduced and full-space operations.
 *
 * @see src/constraints/constraint_transformer.cpp
 * @see src/constraints/builder/builder.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "builder/builder.h"
#include "constraint_map.h"
#include "constraint_set.h"

#include <limits>

namespace fem {
namespace constraint {

/**
 * @class ConstraintTransformer
 * @brief High-level helper that creates and applies constraint maps.
 */
class ConstraintTransformer {
public:
    /**
     * @enum Method
     * @brief Constraint handling backend used by this transformer.
     */
    enum class Method {
        NullSpace, ///< Affine null-space projection: u = u_p + T q.
        Lagrange   ///< Saddle-point system with Lagrange multipliers.
    };

    /**
     * @struct BuildOptions
     * @brief Aggregates options for constraint set assembly and map building.
     */
    struct BuildOptions {
        Method method = Method::NullSpace;   ///< Backend used by the transformer.
        ConstraintSet::Options set;           ///< Options forwarded to `ConstraintSet`.
        ConstraintBuilder::Options builder;   ///< Options forwarded to the builder.
        Precision lagrange_regularization_rel = 1e-6; ///< Relative LM regularization factor (0 disables regularization).
    };

    /**
     * @brief Constructs the transformer and prepares the constraint map.
     *
     * @param eqs Constraint equations to assemble.
     * @param dofs DOF numbering map providing global IDs.
     * @param n_dofs Total number of DOFs.
     */
    ConstraintTransformer(const Equations& eqs,
                          const SystemDofIds& dofs,
                          Index n_dofs);

    /**
     * @brief Constructs the transformer and prepares the constraint map.
     *
     * @param eqs Constraint equations to assemble.
     * @param dofs DOF numbering map providing global IDs.
     * @param n_dofs Total number of DOFs.
     * @param opt Build options for set assembly and map creation.
     */
    ConstraintTransformer(const Equations& eqs,
                          const SystemDofIds& dofs,
                          Index n_dofs,
                          const BuildOptions& opt);

    /// Provides read-only access to the assembled constraint map.
    const ConstraintMap& map() const { return map_; }

    /// Provides read-only access to the assembled constraint set.
    const ConstraintSet& set() const { return set_; }

    /// Provides read-only access to the builder report.
    const ConstraintBuilder::Report& report() const { return report_; }

    /// Returns the active backend.
    Method method() const { return method_; }

    /// Returns a human-readable backend label.
    const char* method_name() const {
        return (method_ == Method::NullSpace) ? "NULLSPACE" : "LAGRANGE";
    }

    /// Indicates whether the constraint system is homogeneous (`d = 0`).
    bool homogeneous() const { return report_.homogeneous; }

    /// Indicates whether the constraint system was deemed feasible.
    bool feasible() const { return report_.feasible; }

    /// Reports the numerical rank determined during factorisation.
    Index rank() const { return report_.rank; }

    /// Reports the number of master DOFs in the reduced system.
    Index n_master() const { return map_.n_master(); }

    /// Assembles the reduced stiffness matrix `A = T^T K T`.
    SparseMatrix assemble_A(const SparseMatrix& K) const;

    /// Assembles the reduced geometric matrix `B = T^T K_g T`.
    SparseMatrix assemble_B(const SparseMatrix& Kg) const;

    /// Assembles the reduced right-hand side `b = T^T (f - K u_p)`.
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const;

    /// Constructs a matrix-free operator for `A = T^T K T`.
    ConstraintMap::OpA opA(const SparseMatrix& K) const;

    /// Constructs a matrix-free operator for `B = T^T K_g T`.
    ConstraintMap::OpB opB(const SparseMatrix& Kg) const;

    /// Applies the transformation `u = T q + u_p`.
    void apply_T(const DynamicVector& q, DynamicVector& u) const;

    /// Applies the adjoint transformation `z = T^T y`.
    void apply_Tt(const DynamicVector& y, DynamicVector& z) const;

    /// Recovers the full solution from reduced coordinates.
    DynamicVector recover_u(const DynamicVector& q) const;

    /// Recovers full-space velocity from reduced coordinates.
    DynamicVector recover_v(const DynamicVector& qdot) const;

    /// Recovers full-space acceleration from reduced coordinates.
    DynamicVector recover_a(const DynamicVector& qddot) const;

    /// Computes reaction forces in the constrained DOFs.
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const;

    /**
     * @brief Estimates Lagrange multipliers for the active constraints.
     *
     * Given a converged static state (u = T q + u_p), this computes a
     * least-squares solution of C^T λ ≈ f - K u, which corresponds to
     * enforcing global equilibrium K u - f + C^T λ = 0 in the presence of
     * eliminated constraints. The returned vector λ has one entry per
     * assembled constraint row (i.e. matches set().m and set().kept_row_ids).
     *
     * @param K Full-space stiffness matrix.
     * @param f Full-space external load vector.
     * @param q Reduced coordinates (solution in master space).
     * @return DynamicVector Lagrange multipliers λ (size = set().m).
     */
    DynamicVector lagrange_multipliers(const SparseMatrix& K,
                                       const DynamicVector& f,
                                       const DynamicVector& q) const;

    /**
     * @brief Builds full-space constraint forces g = C^T λ.
     *
     * @param lambda Lagrange multipliers (size = set().m).
     * @return DynamicVector Full-space constraint force vector.
     */
    DynamicVector constraint_forces(const DynamicVector& lambda) const;

    /**
     * @brief Computes support-only reaction forces from constraint multipliers.
     *
     * Constructs g_supp = C_supp^T λ_supp by accumulating only rows of C that
     * originated from supports. This yields non-zero entries exclusively on DOFs
     * where supports act, even if those DOFs have no element stiffness.
     *
     * @param K Full-space stiffness matrix.
     * @param f Full-space load vector.
     * @param q Reduced coordinates.
     * @return DynamicVector Full-space vector with reactions due to supports only
     *         (sign matches prior convention r = K u - f on supported DOFs).
     */
    DynamicVector support_reactions(const SparseMatrix& K,
                                    const DynamicVector& f,
                                    const DynamicVector& q) const;

    /**
     * @brief Performs diagnostic checks on a static equilibrium solution.
     *
     * @param K Full-space stiffness matrix.
     * @param f Full-space load vector.
     * @param u Full-space displacement vector.
     * @param tol_constraint_rel Relative tolerance for constraint satisfaction.
     * @param tol_reduced_rel Relative tolerance for reduced equilibrium.
     * @param tol_full_rel Relative tolerance for full equilibrium (use `inf` to skip).
     */
    void post_check_static(const SparseMatrix& K,
                           const DynamicVector& f,
                           const DynamicVector& u,
                           Precision tol_constraint_rel = 1e-10,
                           Precision tol_reduced_rel = 1e-8,
                           Precision tol_full_rel = std::numeric_limits<Precision>::infinity()) const;

private:
    void initialize_identity_map(Index n);
    ConstraintBuilder::Report build_lagrange_report() const;
    DynamicVector extract_u_from_solution(const DynamicVector& x) const;
    DynamicVector extract_lambda_from_solution(const DynamicVector& x) const;
    void cache_lagrange_lambda_from_solution(const DynamicVector& x) const;
    DynamicVector solve_multipliers_from_u(const SparseMatrix& K,
                                           const DynamicVector& f,
                                           const DynamicVector& u) const;
    DynamicVector accumulate_support_constraint_forces(const DynamicVector& lambda) const;
    Precision lagrange_regularization_abs(const SparseMatrix& K) const;

    Method method_ = Method::NullSpace; ///< Active constraint backend.
    Precision lagrange_regularization_rel_ = 0; ///< Relative LM regularization (0 disables regularization).
    DynamicVector lagrange_regularization_row_scale_; ///< Per-row scaling for LM regularization (squared row norm of C).
    mutable bool cached_lagrange_lambda_valid_ = false; ///< True if `cached_lagrange_lambda_` matches the most recent LAGRANGE solution.
    mutable DynamicVector cached_lagrange_lambda_; ///< Cached λ extracted from solved LAGRANGE state (size m).
    ConstraintSet set_;               ///< Assembled constraint equations.
    ConstraintMap map_;               ///< Null-space representation.
    ConstraintBuilder::Report report_;///< Diagnostics from the builder.
};

} // namespace constraint
} // namespace fem
