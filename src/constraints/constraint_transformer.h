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
     * @struct BuildOptions
     * @brief Aggregates options for constraint set assembly and map building.
     */
    struct BuildOptions {
        ConstraintSet::Options set;           ///< Options forwarded to `ConstraintSet`.
        ConstraintBuilder::Options builder;   ///< Options forwarded to the builder.
    };

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
                          const BuildOptions& opt = {});

    /// Provides read-only access to the assembled constraint map.
    const ConstraintMap& map() const { return map_; }

    /// Provides read-only access to the assembled constraint set.
    const ConstraintSet& set() const { return set_; }

    /// Provides read-only access to the builder report.
    const ConstraintBuilder::Report& report() const { return report_; }

    /// Indicates whether the constraint system is homogeneous (`d = 0`).
    bool homogeneous() const { return report_.homogeneous; }

    /// Indicates whether the constraint system was deemed feasible.
    bool feasible() const { return report_.feasible; }

    /// Reports the numerical rank determined during factorisation.
    Index rank() const { return report_.rank; }

    /// Reports the number of master DOFs in the reduced system.
    Index n_master() const { return map_.n_master(); }

    /// Assembles the reduced stiffness matrix `A = T^T K T`.
    SparseMatrix assemble_A(const SparseMatrix& K) const { return map_.assemble_A(K); }

    /// Assembles the reduced geometric matrix `B = T^T K_g T`.
    SparseMatrix assemble_B(const SparseMatrix& Kg) const { return map_.assemble_B(Kg); }

    /// Assembles the reduced right-hand side `b = T^T (f - K u_p)`.
    DynamicVector assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
        return map_.assemble_b(K, f);
    }

    /// Constructs a matrix-free operator for `A = T^T K T`.
    ConstraintMap::OpA opA(const SparseMatrix& K) const { return ConstraintMap::OpA(map_, K); }

    /// Constructs a matrix-free operator for `B = T^T K_g T`.
    ConstraintMap::OpB opB(const SparseMatrix& Kg) const { return ConstraintMap::OpB(map_, Kg); }

    /// Applies the transformation `u = T q + u_p`.
    void apply_T(const DynamicVector& q, DynamicVector& u) const { map_.apply_T(q, u); }

    /// Applies the adjoint transformation `z = T^T y`.
    void apply_Tt(const DynamicVector& y, DynamicVector& z) const { map_.apply_Tt(y, z); }

    /// Recovers the full solution from reduced coordinates.
    DynamicVector recover_u(const DynamicVector& q) const { return map_.recover_u(q); }

    /// Computes reaction forces in the constrained DOFs.
    DynamicVector reactions(const SparseMatrix& K, const DynamicVector& f, const DynamicVector& q) const {
        return map_.reactions(K, f, q);
    }

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
    ConstraintSet set_;               ///< Assembled constraint equations.
    ConstraintMap map_;               ///< Null-space representation.
    ConstraintBuilder::Report report_;///< Diagnostics from the builder.
};

} // namespace constraint
} // namespace fem
