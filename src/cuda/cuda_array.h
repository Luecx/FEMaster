/**
 * @file cuda_array.h
 * @brief Declares a simple RAII wrapper for device memory allocations.
 *
 * @see src/cuda/cuda_array.cu
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#ifdef SUPPORT_GPU

#include "cuda.h"
#include "assert_cuda.h"
#include "../core/assert.h"

#include <type_traits>
#include <vector>

namespace fem {
namespace cuda {

template<typename T>
class CudaArray {
public:
    explicit CudaArray(Index size);
    ~CudaArray();

    operator T*();

    Index size() const;

    void copy(const T* source);

    template<typename SourceType>
    void upload(const SourceType* source);

    template<typename DestType>
    void download(DestType* destination);

    void clear();

    static Index estimate_mem(Index count);

protected:
    void allocate();
    void release();

    T* data = nullptr;
    Index element_count = 0;
};
} // namespace cuda
} // namespace fem

#include "cuda_array.tpp"

#endif
