/******************************************************************************
* @file eigen_internal.h
 * @brief Internal entry points for the eigen solver implementations.
 *
 * @details These functions implement the individual solver flavors and are
 *          called from the public dispatcher in `eigen.cpp`. They assume
 *          all inputs have already been validated and compressed.
 *
 *          Implementations live in separate translation units for faster
 *          incremental builds and clearer code ownership per mode.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 ******************************************************************************/
#pragma once
#include "eigen_common.h"

namespace fem::solver::detail {

/** Simple: A x = λ x (Regular) */
std::vector<EigenValueVectorPair>
eigs_simple_regular(const SparseMatrix& A, int k, const EigenOpts& opts);

/** Simple: A x = λ x (ShiftInvert via (A - σI)^{-1}) */
std::vector<EigenValueVectorPair>
eigs_simple_shiftinvert(const SparseMatrix& A, int k, const EigenOpts& opts);

/** Generalized: A x = λ B x (ShiftInvert) */
std::vector<EigenValueVectorPair>
eigs_general_shiftinvert(const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts);

/** Generalized: A x = λ B x (Buckling) */
std::vector<EigenValueVectorPair>
eigs_general_buckling(const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts);

/** Generalized: A x = λ B x (Cayley) */
std::vector<EigenValueVectorPair>
eigs_general_cayley(const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts);

} // namespace fem::solver::detail
