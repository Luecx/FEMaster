/**
 * @file constraint_transformer.h
 * @brief Declares the common interface for assembling and enforcing linear constraints.
 *
 * `ConstraintTransformer` converts symbolic constraint equations into the
 * algebraic system
 *
 *     C u = d
 *
 * and applies one of three enforcement strategies:
 *
 * - affine null-space reduction,
 * - direct elimination,
 * - Lagrange multipliers.
 *
 * For reduced methods, admissible full-system vectors are represented as
 *
 *     u = u_p + T q,
 *
 * where `u_p` is a particular solution of the inhomogeneous constraints, `T`
 * spans the homogeneous null space of `C` and `q` contains the independent
 * reduced coordinates.
 *
 * For the Lagrange formulation, the original displacement unknowns are retained
 * and augmented by one multiplier per assembled constraint equation.
 *
 * Besides system assembly and coordinate recovery, the class provides reaction
 * reconstruction and diagnostic post-checks for constraint satisfaction and
 * equilibrium.
 *
 * @see constraint_transformer_core.cpp
 * @see constraint_transformer_lagrange.cpp
 * @see constraint_transformer_postcheck.cpp
 * @see constraint_system.h
 * @see constraint_map.h
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "constraint_build_report.h"
#include "constraint_map.h"
#include "constraint_system.h"
#include "null_space.h"

#include <limits>

namespace fem {
namespace constraint {

/**
 * @brief Assembles constraints and exposes a method-independent solver interface.
 *
 * The transformer owns the assembled constraint system and the representation
 * required by the selected enforcement method.
 *
 * Null-space reduction and elimination generate a `ConstraintMap` describing
 *
 *     u = u_p + T q.
 *
 * The Lagrange method instead constructs the augmented saddle-point system
 *
 *     [ K   C^T ] [ u      ] = [ f ]
 *     [ C    0  ] [ lambda ]   [ d ],
 *
 * with optional row scaling and diagonal regularization of the multiplier
 * block.
 *
 * Public recovery functions distinguish affine displacement quantities from
 * homogeneous increments. Displacements include the particular solution
 * `u_p`, whereas increments, velocities and accelerations are transformed only
 * through the homogeneous matrix `T`.
 */
class ConstraintTransformer {


public:
    /**
     * @brief Available strategies for enforcing the assembled constraint equations.
     */
    enum class Method {
        // Construct an affine null-space basis from C and solve only for the
        // independent reduced coordinates q.
        NullSpace,

        // Retain the full displacement vector and augment the algebraic system
        // with one Lagrange multiplier per constraint equation.
        Lagrange,

        // Select one dependent DOF per independent equation and recursively
        // substitute it through the remaining master DOFs.
        Elimination
    };

    /**
     * @brief Configures constraint assembly and the selected enforcement method.
     */
    struct Options {
        // Constraint enforcement strategy used by the transformer.
        Method method{};

        // Options controlling conversion of symbolic equations into the sparse
        // global system C u = d.
        ConstraintOptions constraints{};

        // Numerical settings used only when Method::NullSpace is selected.
        NullSpaceOptions null_space{};

        // Relative regularization applied to the diagonal of the Lagrange
        // multiplier block. The effective value is scaled with a characteristic
        // stiffness extracted from the supplied system matrix.
        Precision lagrange_regularization{Precision(1e-10)};
    };

private:
    // Selected enforcement strategy.
    Method method_{};

    // Non-negative relative regularization configured for the Lagrange
    // multiplier block.
    Precision lagrange_regularization_{Precision(0)};

    // Inverse row norms used to scale the Lagrange constraint equations.
    DynamicVector lagrange_row_scale_{};

    // Assembled global constraint matrix, right-hand side and source metadata.
    ConstraintSystem system_{};

    // Affine reduced-coordinate representation used by null-space reduction and
    // direct elimination.
    ConstraintMap map_{};

    // Numerical diagnostics describing rank, redundancy, homogeneity and
    // feasibility of the constructed representation.
    ConstraintBuildReport report_{};
public:

    // Construct a transformer using default options. The default-initialized
    // method stored in Options determines the enforcement strategy.
    ConstraintTransformer(const Equations&     equations,
                          const SystemDofIds& dofs,
                          Index                n_dofs);

    // Assemble the symbolic equations and initialize the explicitly selected
    // constraint representation.
    ConstraintTransformer(const Equations&     equations,
                          const SystemDofIds& dofs,
                          Index                n_dofs,
                          const Options&       options);

    // Return the selected constraint enforcement method.
    Method method() const { return method_; }

    // Return a stable uppercase name for the selected method, suitable for
    // logging and diagnostic output.
    const char* method_name() const;

    // Return the assembled sparse constraint system C u = d.
    const ConstraintSystem& system() const { return system_; }

    // Return the diagnostics generated while constructing the selected
    // constraint representation.
    const ConstraintBuildReport& report() const { return report_; }

    // Return the affine reduced-coordinate map. This function is valid only for
    // null-space reduction and elimination.
    const ConstraintMap& reduced_map() const;

    // Return true when the selected method owns an affine reduced-coordinate
    // map. The Lagrange formulation does not construct one.
    bool has_reduced_map() const { return method_ != Method::Lagrange; }

    // Return true when the assembled right-hand side d is zero.
    bool homogeneous() const { return report_.homogeneous; }

    // Return true when the generated representation satisfies the assembled
    // constraints within the applicable feasibility tolerance.
    bool feasible() const { return report_.feasible; }

    // Return true when the numerical rank has been determined explicitly.
    bool rank_known() const { return report_.rank_known; }

    // Return the numerical rank reported for the assembled constraint system.
    Index rank() const { return report_.rank; }

    // Return the number of algebraic unknowns in the transformed solver system.
    //
    // Reduced methods return the number of master coordinates, while the
    // Lagrange method returns the displacement DOFs plus the multipliers.
    Index unknowns() const;

    // Return the number of independent reduced displacement coordinates.
    // This operation is unavailable for the Lagrange formulation.
    Index independent_dofs() const;

    // Assemble the primary matrix used by the selected solver formulation.
    //
    // Reduced methods return T^T K T. The Lagrange method returns the augmented
    // saddle-point matrix containing K and the scaled constraint blocks.
    SparseMatrix assemble_system_matrix(const SparseMatrix& matrix) const;

    // Reduce a secondary displacement-space matrix through T^T A T. This is
    // intended for matrices such as mass or geometric stiffness and is
    // unavailable for the Lagrange formulation.
    SparseMatrix reduce_secondary_matrix(const SparseMatrix& matrix) const;

    // Assemble the right-hand side matching the selected primary system.
    //
    // Reduced methods return T^T (f - K u_p). The Lagrange method appends the
    // scaled constraint right-hand side to the original load vector.
    DynamicVector assemble_system_rhs(const SparseMatrix&  matrix,
                                      const DynamicVector& rhs) const;

    // Apply T^T to a full-system vector and write the result into the supplied
    // reduced vector. This operation is available only for reduced methods.
    void project_vector(const DynamicVector& full,
                        DynamicVector&       reduced) const;

    // Allocate and return the reduced projection T^T full.
    DynamicVector project_vector(const DynamicVector& full) const;

    // Recover the physical displacement vector represented by a solver
    // solution. Reduced methods evaluate u_p + T q; the Lagrange method extracts
    // the displacement block from the augmented solution.
    DynamicVector recover_displacement(const DynamicVector& solution) const;

    // Recover a homogeneous displacement increment. The particular solution is
    // intentionally excluded because increments represent directions rather
    // than absolute affine configurations.
    DynamicVector recover_increment(const DynamicVector& increment) const;

    // Recover a full-system velocity from reduced or augmented coordinates.
    // Velocities use the homogeneous transformation and therefore follow the
    // same path as displacement increments.
    DynamicVector recover_velocity(const DynamicVector& velocity) const;

    // Recover a full-system acceleration from reduced or augmented coordinates.
    // Accelerations also exclude the particular displacement contribution.
    DynamicVector recover_acceleration(const DynamicVector& acceleration) const;

    // Compute the complete full-system residual
    //
    //     r = K u - f
    //
    // after recovering the physical displacement vector from the supplied
    // solver solution.
    DynamicVector reactions(const SparseMatrix&  matrix,
                            const DynamicVector& rhs,
                            const DynamicVector& solution) const;

    // Recover only the force contribution associated with equations originating
    // from supports. For Lagrange systems, the stored multiplier block is used
    // directly. Reduced methods reconstruct compatible multipliers from the
    // equilibrium residual.
    DynamicVector support_reactions(const SparseMatrix&  matrix,
                                    const DynamicVector& rhs,
                                    const DynamicVector& solution) const;

    // Evaluate static-solution diagnostics for the selected formulation.
    //
    // The checks include constraint satisfaction and equilibrium in the
    // relevant coordinate space. Reduced methods additionally verify affine
    // consistency of u_p and sample the homogeneous invariant C T = 0.
    //
    // The full residual check is disabled by default through an infinite
    // tolerance because constrained systems generally retain reaction forces in
    // the full-space residual.
    void post_check_static(
        const SparseMatrix&  matrix,
        const DynamicVector& rhs,
        const DynamicVector& solution,
        Precision            constraint_tolerance = Precision(1e-10),
        Precision            reduced_tolerance    = Precision(1e-8),
        Precision            full_tolerance       = std::numeric_limits<Precision>::infinity()
    ) const;

private:
    // Initialize row scaling and diagnostic metadata required by the Lagrange
    // multiplier formulation.
    void initialize_lagrange();

    // Assemble the augmented Lagrange matrix from the supplied displacement-
    // space matrix and the scaled constraint matrix.
    SparseMatrix assemble_lagrange_matrix(const SparseMatrix& matrix) const;

    // Assemble the augmented Lagrange right-hand side containing the physical
    // load vector and scaled prescribed constraint values.
    DynamicVector assemble_lagrange_rhs(const DynamicVector& rhs) const;

    // Extract the physical displacement block from an augmented Lagrange
    // solution vector.
    DynamicVector extract_lagrange_displacement(const DynamicVector& solution) const;

    // Extract the multiplier block from an augmented Lagrange solution vector.
    DynamicVector extract_lagrange_multipliers(const DynamicVector& solution) const;

    // Reconstruct constraint multipliers from a recovered displacement and the
    // full equilibrium residual. This path is used for reduced methods when
    // source-specific constraint forces are requested.
    DynamicVector solve_multipliers(const SparseMatrix&  matrix,
                                    const DynamicVector& rhs,
                                    const DynamicVector& displacement) const;

    // Convert multiplier values into a full-system force vector while retaining
    // only rows originating from support equations.
    DynamicVector support_constraint_forces(const DynamicVector& multipliers,
                                            bool                 scaled_rows) const;

    // Apply the precomputed Lagrange row scaling to a vector defined over the
    // constraint equations.
    DynamicVector scale_lagrange_rows(const DynamicVector& values) const;

    // Compute the absolute multiplier-block regularization from the configured
    // relative factor and a characteristic scale of the supplied system matrix.
    Precision lagrange_regularization(const SparseMatrix& matrix) const;
};

} // namespace constraint
} // namespace fem

