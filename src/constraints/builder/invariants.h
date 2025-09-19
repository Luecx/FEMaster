/******************************************************************************
* @file invariants.h
 * @brief Optional debug-time invariant checks for T and fast apply methods.
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"
#include "../constraint_map.h"

namespace fem { namespace constraint {

void check_invariants(const ConstraintMap& M, const SparseMatrix& C);

}} // namespace
