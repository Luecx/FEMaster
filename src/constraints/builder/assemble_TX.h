/**
 * @file assemble_TX.h
 * @brief Declares helpers that assemble the `T` and `X` matrices.
 *
 * The assembly combines reduced column data and partitioning information to
 * build the null-space transformation matrices used by `ConstraintMap`.
 *
 * @see src/constraints/builder/assemble_TX.cpp
 * @see src/constraints/builder/build_x.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"
#include "build_x.h"

#include <vector>

namespace fem {
namespace constraint {

/**
 * @struct AssembleInput
 * @brief Holds the data required to assemble `T` and `X`.
 */
struct AssembleInput {
    int n = 0;                       ///< Total number of DOFs.
    int r = 0;                       ///< Numerical rank.
    std::vector<int> slaves_loc;     ///< Local slave indices.
    std::vector<int> masters_loc;    ///< Local master indices.
    std::vector<int> used;           ///< Mapping from local to global DOFs.
    std::vector<Index> masters_glob; ///< Final master ordering in global indices.
};

/**
 * @struct AssembleOutput
 * @brief Contains the assembled `T` and `X` matrices.
 */
struct AssembleOutput {
    SparseMatrix T; ///< Transformation matrix mapping reduced to full DOFs.
    SparseMatrix X; ///< Constraint matrix expressing slaves in terms of masters.
};

/**
 * @brief Builds the `T` and `X` matrices from partitioning and column data.
 *
 * @param in Assembly input data.
 * @param cols Reduced column data for `X`.
 * @return AssembleOutput Populated matrices.
 */
AssembleOutput assemble_T_and_X(const AssembleInput& in, const XCols& cols);

} // namespace constraint
} // namespace fem
