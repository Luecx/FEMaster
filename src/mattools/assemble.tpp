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

    Eigen::setNbThreads(1);  // Ensure Eigen only uses one thread per operation

    constexpr int MAX_LOCAL_SIZE = 128;
    constexpr int BASE_BATCH_SIZE = 1024 * 1024 * 16;

    SparseMatrix matrix{indices.maxCoeff() + 1, indices.maxCoeff() + 1};

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = global_config.max_threads;
#endif

    int adjusted_batch_size = BASE_BATCH_SIZE / num_threads;

    // Create a vector of vectors for triplets, one per thread
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_per_thread(num_threads);
    for (auto &triplets : triplets_per_thread) {
        triplets.reserve(adjusted_batch_size);
    }

    alignas(64) double local_storage[MAX_LOCAL_SIZE * MAX_LOCAL_SIZE]{};

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];
        if (element == nullptr) continue;

        // Compute the thread ID
#ifdef _OPENMP
        int thread_id = omp_get_thread_num();
#else
        int thread_id = 0;
#endif

        auto& local_triplets = triplets_per_thread[thread_id];

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

                        // Add the local matrix contribution to the thread-local triplets
                        local_triplets.emplace_back(global_i, global_j, local_matrix(i * el_ndofs + idof, j * el_ndofs + jdof));
                    }
                }
            }
        }

        // Flush triplets to prevent excessive memory use
        if (local_triplets.size() > adjusted_batch_size * 9 / 10) {
            #pragma omp critical
            {
                matrix.insertFromTriplets(local_triplets.begin(), local_triplets.end());
                local_triplets.clear();
            }
        }
    }

    // Final flush of all triplets from all threads
    for (auto& local_triplets : triplets_per_thread) {
        matrix.insertFromTriplets(local_triplets.begin(), local_triplets.end());
    }

    return matrix;
}

} } // namespace fem::mattools
