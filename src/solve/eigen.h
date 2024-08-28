#pragma once

#include "../cuda/cuda_csr.h"
#include "../cuda/cuda_array.h"
#include "../cuda/cuda_vec.h"
#include "../core/core.h"
#include "device.h"
#include "method.h"

#include <string>
#include <utility>  // For std::pair

namespace solver {

/******************************************************************************
 * @brief Extracts eigenvalues (and optionally eigenvectors) from a matrix using a method on the specified device.
 *
 * @param device The computational device to use (CPU or GPU).
 * @param mat The sparse matrix for which eigenvalues are to be computed.
 * @param num_eigenvalues The number of eigenvalues (and eigenvectors) to compute.
 * @param return_vectors If true, the function returns both eigenvalues and eigenvectors.
 *                       If false, only the eigenvalues are returned.
 * @return std::pair<DynamicVector, DynamicMatrix> A pair containing the eigenvalues and (optionally) eigenvectors.
 *         If `return_vectors` is false, the second element of the pair will be an empty matrix.
 ******************************************************************************/
std::pair<DynamicVector, DynamicMatrix> compute_eigenvalues(SolverDevice device,
                                                            SparseMatrix& mat,
                                                            int num_eigenvalues,
                                                            bool return_vectors = false);

}  // namespace solver
