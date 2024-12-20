/******************************************************************************
 * @file assemble.h
 * @brief Provides assembly functions for system matrices in FEM using
 * user-defined lambda functions and element pointers.
 *
 * The user provides a lambda that handles the local matrix contributions for
 * each element, including the local matrix computations. The lambda function
 * is called for each element pointer, and the global matrix is assembled accordingly.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include <Eigen/Sparse>
#include <vector>
#include <memory>  // for std::shared_ptr

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Assembles a global sparse matrix in FEM using a user-provided lambda
 * function that generates local element contributions for each element.
 *
 * This function iterates over a given vector of element pointers and calls the
 * lambda function for each element to obtain the local element matrix. The local
 * matrix is then assembled into the global sparse matrix using the provided index
 * matrix that maps local element DOFs to global DOFs.
 *
 * The lambda must have the following signature:
 *
 * `MapMatrix lambda(ElementPtr element, Precision* local_storage)`
 *
 * The lambda is responsible for generating the local element matrix based on the
 * provided `element`, a preallocated local storage array, and the coordinates of
 * the nodes.
 *
 * @tparam Lambda A callable object (e.g., a lambda function) that generates local
 * matrix contributions for each element.
 *
 * @param elements A vector of element pointers that will be processed.
 * @param indices An SystemDofIds mapping element node and DOF indices to global DOF indices.
 * @param lambda The lambda function responsible for producing the local element matrix.
 * @return SparseMatrix The assembled global sparse matrix.
 ******************************************************************************/
template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements, const SystemDofIds &indices, Lambda&& lambda);

} } // namespace fem::mattools

#include "assemble.tpp"