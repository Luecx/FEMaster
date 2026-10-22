/**
 * @file build_sparsity_pattern.h
 * @brief Declares connectivity-based sparse matrix pattern construction.
 *
 * The matrix tools layer derives admissible global matrix positions from element
 * connectivity and active nodal DOF numbering. The result is a compressed sparse
 * matrix with zero-valued coefficients, ready for numerical assembly.
 *
 * Element matrix evaluation, material-state updates and nodal force accumulation
 * remain responsibilities of assemble_matrix and its callbacks. Pattern caching
 * and reuse belong to the caller; no model or material state is changed here.
 *
 * @see build_sparsity_pattern
 * @see assemble_matrix
 *
 * @author Finn Eggers
 * @date 22.10.2026
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"

namespace fem { namespace mattools {

// Build a conservative compressed-column pattern from connected nodes and all
// their active DOFs. Return zero-valued coefficients without evaluating local matrices.
// Active global indices must be unique; negative indices are inactive.
SparseMatrix build_sparsity_pattern(const std::vector<model::ElementPtr>& elements,
                                    const SystemDofIds& indices);

} } // namespace fem::mattools
