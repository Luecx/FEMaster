/******************************************************************************
 * @file assemble.cpp
 * @brief Provides the implementation of the assemble_matrix function for system
 * matrices in FEM, using user-defined lambda functions.
 *
 * The assembly process handles the global sparse matrix construction by iterating
 * over a vector of element pointers and calling the lambda for local matrix
 * contributions, which are then assembled into the global matrix.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 ******************************************************************************/

namespace fem { namespace mattools {

template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements, const SystemDofIds &indices, Lambda&& lambda) {
    constexpr int BATCH_SIZE = 1024 * 1024 * 16;
    constexpr int MAX_LOCAL_SIZE = 128;

    SparseMatrix matrix{indices.maxCoeff() + 1, indices.maxCoeff() + 1};
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(BATCH_SIZE);

    alignas(64) double local_storage[MAX_LOCAL_SIZE * MAX_LOCAL_SIZE]{};

    for (const auto& element : elements) {
        if (element == nullptr) continue;

        // Call the lambda to get the local element matrix (K) for the current element
        auto local_matrix = lambda(element, local_storage);

        // Element-specific properties
        auto el_nodes = element->n_nodes();
        auto el_ndofs = local_matrix.rows() / el_nodes;

        for (int i = 0; i < el_nodes; i++) {
            for (int j = 0; j < el_nodes; j++) {
                for (int idof = 0; idof < el_ndofs; idof++) {
                    for (int jdof = 0; jdof < el_ndofs; jdof++) {
                        int global_i = indices(element->nodes()[i], idof);
                        int global_j = indices(element->nodes()[j], jdof);

                        // Add the local matrix contribution to the global matrix using triplets
                        triplets.emplace_back(global_i, global_j, local_matrix(i * el_ndofs + idof, j * el_ndofs + jdof));
                    }
                }
            }
        }

        // Periodically flush triplets to the sparse matrix to prevent excessive memory use
        if (triplets.size() > BATCH_SIZE * 9 / 10) {
            matrix.insertFromTriplets(triplets.begin(), triplets.end());
            triplets.clear();
        }
    }

    // Final flush of the remaining triplets
    matrix.insertFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}

} } // namespace fem::mattools
