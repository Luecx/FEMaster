/**
 * @file constraint_transformer.h
 * @brief Declares the common interface for constraint transformation methods.
 *
 * @see src/constraints/transformer/constraint_transformer_core.cpp
 * @see src/constraints/transformer/constraint_transformer_lagrange.cpp
 * @see src/constraints/transformer/constraint_transformer_postcheck.cpp
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
 * @brief Builds and applies the selected constraint transformation method.
 */
class ConstraintTransformer {
public:
    /**
     * @brief Available methods for enforcing constraint equations.
     */
    enum class Method {
        NullSpace,  ///< Null-space reduction using `u = u_p + T q`.
        Lagrange,   ///< Lagrange multiplier formulation.
        Elimination ///< Direct elimination of dependent DOFs.
    };

    /**
     * @brief Configuration used to construct the constraint transformer.
     */
    struct Options {
        // General constraint settings
        Method            method{};      ///< Constraint enforcement method.
        ConstraintOptions constraints{}; ///< Constraint matrix assembly options.

        // Method-specific settings
        NullSpaceOptions null_space{}; ///< Null-space construction options.

        Precision lagrange_regularization{Precision(1e-10)}; ///< Relative regularization applied to the Lagrange block.
    };

    // Construction
    ConstraintTransformer(const Equations&     equations,
                          const SystemDofIds& dofs,
                          Index                n_dofs);

    ConstraintTransformer(const Equations&     equations,
                          const SystemDofIds& dofs,
                          Index                n_dofs,
                          const Options&       options);

    // Method information
    Method      method() const { return method_; }
    const char* method_name() const;

    // Constraint data
    const ConstraintSystem&      system() const { return system_; }
    const ConstraintBuildReport& report() const { return report_; }
    const ConstraintMap&         reduced_map() const;

    // Constraint properties
    bool has_reduced_map() const { return method_ != Method::Lagrange; }
    bool homogeneous()     const { return report_.homogeneous; }
    bool feasible()        const { return report_.feasible; }
    bool rank_known()      const { return report_.rank_known; }

    // System dimensions
    Index rank() const { return report_.rank; }
    Index unknowns() const;
    Index independent_dofs() const;

    // System assembly
    SparseMatrix  assemble_system_matrix(const SparseMatrix& matrix) const;
    SparseMatrix  reduce_secondary_matrix(const SparseMatrix& matrix) const;
    DynamicVector assemble_system_rhs(const SparseMatrix&  matrix,
                                      const DynamicVector& rhs) const;

    // Coordinate transformations
    void project_vector(const DynamicVector& full,
                        DynamicVector&       reduced) const;

    DynamicVector project_vector      (const DynamicVector& full) const;
    DynamicVector recover_displacement(const DynamicVector& solution) const;
    DynamicVector recover_increment   (const DynamicVector& increment) const;
    DynamicVector recover_velocity    (const DynamicVector& velocity) const;
    DynamicVector recover_acceleration(const DynamicVector& acceleration) const;

    // Reaction forces
    DynamicVector reactions(const SparseMatrix&  matrix,
                            const DynamicVector& rhs,
                            const DynamicVector& solution) const;

    DynamicVector support_reactions(const SparseMatrix&  matrix,
                                    const DynamicVector& rhs,
                                    const DynamicVector& solution) const;

    // Diagnostic checks
    void post_check_static(
        const SparseMatrix&  matrix,
        const DynamicVector& rhs,
        const DynamicVector& solution,
        Precision            constraint_tolerance = Precision(1e-10),
        Precision            reduced_tolerance    = Precision(1e-8),
        Precision            full_tolerance       = std::numeric_limits<Precision>::infinity()
    ) const;

private:
    // Lagrange initialization and assembly
    void initialize_lagrange();

    SparseMatrix  assemble_lagrange_matrix(const SparseMatrix& matrix) const;
    DynamicVector assemble_lagrange_rhs   (const DynamicVector& rhs) const;

    // Lagrange solution extraction
    DynamicVector extract_lagrange_displacement(const DynamicVector& solution) const;
    DynamicVector extract_lagrange_multipliers (const DynamicVector& solution) const;

    // Reaction recovery
    DynamicVector solve_multipliers(const SparseMatrix&  matrix,
                                    const DynamicVector& rhs,
                                    const DynamicVector& displacement) const;

    DynamicVector support_constraint_forces(const DynamicVector& multipliers,
                                            bool                 scaled_rows) const;

    // Lagrange scaling and regularization
    DynamicVector scale_lagrange_rows(const DynamicVector& values) const;
    Precision     lagrange_regularization(const SparseMatrix& matrix) const;

    // Selected method and settings
    Method    method_{};
    Precision lagrange_regularization_{Precision(0)};

    // Assembled constraint data
    DynamicVector         lagrange_row_scale_{};
    ConstraintSystem      system_{};
    ConstraintMap         map_{};
    ConstraintBuildReport report_{};
};

} // namespace constraint
} // namespace fem

