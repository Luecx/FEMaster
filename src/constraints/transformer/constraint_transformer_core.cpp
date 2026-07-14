/**
 * @file constraint_transformer_core.cpp
 * @brief Implements the common constraint transformation operations.
 *
 * @see src/constraints/transformer/constraint_transformer.h
 * @see src/constraints/transformer/constraint_transformer_lagrange.cpp
 * @see src/constraints/transformer/constraint_transformer_postcheck.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../../core/logging.h"
#include "elimination.h"

#include <algorithm>
#include <tuple>

namespace fem {
namespace constraint {

ConstraintTransformer::ConstraintTransformer(const Equations&    equations,
                                             const SystemDofIds& dofs,
                                             Index                n_dofs)
    : ConstraintTransformer(equations, dofs, n_dofs, Options{}) {}

ConstraintTransformer::ConstraintTransformer(const Equations&    equations,
                                             const SystemDofIds& dofs,
                                             Index                n_dofs,
                                             const Options&       options)
    : method_                  (options.method),
      lagrange_regularization_ (std::max(Precision(0), options.lagrange_regularization)),
      system_                  (assemble_constraint_system(
          equations,
          dofs,
          n_dofs,
          options.constraints
      )) {
    // Build the selected constraint representation
    if (method_ == Method::NullSpace) {
        std::tie(map_, report_) = build_null_space(system_, options.null_space);
    } else if (method_ == Method::Elimination) {
        std::tie(map_, report_) = build_elimination(system_);
    } else {
        initialize_lagrange();
    }
}

const char* ConstraintTransformer::method_name() const {
    switch (method_) {
        case Method::NullSpace:   return "NULLSPACE";
        case Method::Lagrange:    return "LAGRANGE";
        case Method::Elimination: return "ELIMINATION";
    }

    return "UNKNOWN";
}

const ConstraintMap& ConstraintTransformer::reduced_map() const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] reduced map requested for LAGRANGE"
    );

    return map_;
}

Index ConstraintTransformer::unknowns() const {
    return method_ == Method::Lagrange
        ? system_.dofs + system_.equations
        : map_.n_master();
}

Index ConstraintTransformer::independent_dofs() const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] independent DOFs requested for LAGRANGE"
    );

    return map_.n_master();
}

SparseMatrix ConstraintTransformer::assemble_system_matrix(const SparseMatrix& matrix) const {
    return method_ == Method::Lagrange
        ? assemble_lagrange_matrix(matrix)
        : map_.reduce_matrix(matrix);
}

SparseMatrix ConstraintTransformer::reduce_secondary_matrix(const SparseMatrix& matrix) const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] secondary matrix reduction is unavailable for LAGRANGE"
    );

    return map_.reduce_matrix(matrix);
}

DynamicVector ConstraintTransformer::assemble_system_rhs(const SparseMatrix&  matrix,
                                                         const DynamicVector& rhs) const {
    return method_ == Method::Lagrange
        ? assemble_lagrange_rhs(rhs)
        : map_.reduce_rhs(matrix, rhs);
}

void ConstraintTransformer::project_vector(const DynamicVector& full,
                                           DynamicVector&       reduced) const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] vector projection is unavailable for LAGRANGE"
    );

    map_.project(full, reduced);
}

DynamicVector ConstraintTransformer::project_vector(const DynamicVector& full) const {
    DynamicVector reduced{};
    project_vector(full, reduced);
    return reduced;
}

DynamicVector ConstraintTransformer::recover_displacement(const DynamicVector& solution) const {
    return method_ == Method::Lagrange
        ? extract_lagrange_displacement(solution)
        : map_.recover(solution);
}

DynamicVector ConstraintTransformer::recover_increment(const DynamicVector& increment) const {
    return method_ == Method::Lagrange
        ? extract_lagrange_displacement(increment)
        : map_.T * increment;
}

DynamicVector ConstraintTransformer::recover_velocity(const DynamicVector& velocity) const {
    return recover_increment(velocity);
}

DynamicVector ConstraintTransformer::recover_acceleration(const DynamicVector& acceleration) const {
    return recover_increment(acceleration);
}

DynamicVector ConstraintTransformer::reactions(const SparseMatrix&  matrix,
                                               const DynamicVector& rhs,
                                               const DynamicVector& solution) const {
    return matrix * recover_displacement(solution) - rhs;
}

} // namespace constraint
} // namespace fem

