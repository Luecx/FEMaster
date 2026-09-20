/**
 * @file parallel.h
 * @brief Provides indexed OpenMP execution with safe C++ exception propagation.
 *
 * FEMaster uses independent element and nodal loops in assembly, step
 * initialization and result recovery. This header centralizes the repeated
 * OpenMP scheduling and exception handling without owning the numerical
 * operation, its output fields or its worker-local accumulators.
 *
 * The caller may request a minimum scheduling chunk size to avoid starting
 * many workers for small jobs. With the default minimum of zero, the helper
 * aims for approximately sixteen scheduling chunks per worker. This is
 * distinct from the triplet buffer size used during sparse matrix assembly.
 *
 * @see fem::mattools::assemble_matrix
 */

#pragma once

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>

#include "config.h"

#include <Eigen/Core>

#ifdef USE_MKL
    #include <mkl_service.h>
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem::parallel {

/**
 * Restricts Eigen and MKL to one internal worker for an explicit loop and
 * restores FEMaster's configured thread count on every exit, including throws.
 * Instantiate on the calling thread, outside any OpenMP worker region.
 */
class ScopedLinearAlgebraThreads {
public:
    ScopedLinearAlgebraThreads() { limit(); }

    ScopedLinearAlgebraThreads(const ScopedLinearAlgebraThreads&) = delete;
    ScopedLinearAlgebraThreads& operator=(const ScopedLinearAlgebraThreads&) = delete;

    ~ScopedLinearAlgebraThreads() {
        Eigen::setNbThreads(global_config.max_threads);
#ifdef USE_MKL
        mkl_set_num_threads(global_config.max_threads);
#endif
    }

    // The indexed helper restores the configured counts when it returns.
    // Call this again if the caller has further explicit OpenMP work afterwards.
    void limit() const {
        Eigen::setNbThreads(1);
#ifdef USE_MKL
        mkl_set_num_threads(1);
#endif
    }
};

/**
 * Returns the useful worker count for an independent indexed loop.
 *
 * A nonzero minimum batch size limits the worker count so that each worker
 * receives at least that many iterations. A value of zero imposes no minimum
 * beyond one iteration per worker. Calls inside an existing OpenMP region
 * execute serially to avoid allocating worker-local storage for nested teams.
 *
 * @param count Number of iterations.
 * @param max_threads Upper bound on the number of workers.
 * @param min_batch_size Optional minimum iterations per worker; zero by default.
 * @return Number of workers to allocate, always at least one.
 */
template<typename Index>
int worker_count(Index count, int max_threads, Index min_batch_size = 0) {
#ifdef _OPENMP
    if (count <= 0 || max_threads <= 1 || omp_in_parallel()) {
        return 1;
    }

    // Limit workers by the number of iterations available per worker
    const Index minimum = std::max<Index>(Index(1), min_batch_size);
    const Index useful  = count / minimum;
    return static_cast<int>(std::max<Index>(
        Index(1),
        std::min<Index>(useful, static_cast<Index>(max_threads))
    ));
#else
    (void) count;
    (void) max_threads;
    (void) min_batch_size;
    return 1;
#endif
}

/**
 * Executes independent iterations and propagates worker exceptions safely.
 *
 * The chunk size targets approximately sixteen scheduling chunks per worker
 * while respecting the optional minimum batch size. The final chunk may
 * contain fewer iterations. The callback receives the iteration index and
 * worker index; it must write to disjoint output or use worker-local storage.
 *
 * Exceptions are captured within the worker and rethrown only after all
 * workers leave the OpenMP region. Other iterations may finish after a worker
 * fails; callers remain responsible for cleaning up partially updated state.
 *
 * @param count Number of iterations.
 * @param max_threads Maximum number of workers to use.
 * @param function Callback receiving the iteration index and worker index.
 * @param min_batch_size Optional minimum scheduling chunk size; zero by default.
 */
template<typename Index, typename Function>
void for_index(Index count, int max_threads, Function&& function, Index min_batch_size = 0) {
    if (count <= 0) {
        return;
    }

    // Applies to both the serial and OpenMP paths; RAII also restores the
    // configured Eigen/MKL counts before propagating a callback exception.
    const ScopedLinearAlgebraThreads threading;

#ifdef _OPENMP
    const int threads = worker_count(count, max_threads, min_batch_size);

    if (threads > 1) {
        // Aim for sixteen chunks per worker without shrinking below the minimum
        const Index minimum   = std::max<Index>(Index(1), min_batch_size);
        const Index batch_size = std::max<Index>(
            minimum,
            (count - Index(1)) / static_cast<Index>(threads) / Index(16) + Index(1)
        );

        std::atomic_bool   failed{false};
        std::exception_ptr failure = nullptr;
        std::mutex         failure_mutex;

        // Keep all C++ exceptions inside workers until the OpenMP barrier
        #pragma omp parallel for schedule(static, batch_size) num_threads(threads)
        for (Index i = 0; i < count; ++i) {
            if (failed.load(std::memory_order_relaxed)) {
                continue;
            }

            try {
                function(i, omp_get_thread_num());
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) {
                    failure = std::current_exception();
                }
                failed.store(true, std::memory_order_relaxed);
            }
        }

        if (failure) {
            std::rethrow_exception(failure);
        }

        return;
    }
#else
    (void) max_threads;
    (void) min_batch_size;
#endif

    // Avoid worker-team overhead for small loops and builds without OpenMP
    for (Index i = 0; i < count; ++i) {
        function(i, 0);
    }
}

} // namespace fem::parallel
