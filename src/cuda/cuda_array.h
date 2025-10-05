/******************************************************************************
 * @file cuda_array.h
 * @brief Declares a simple RAII wrapper for device memory allocations.
 *
 * @see src/cuda/cuda_array.cu
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#ifdef SUPPORT_GPU

#include "cuda.h"

#include <type_traits>
#include <vector>

namespace fem {
namespace cuda {

template<typename T>
class CudaArray {
public:
    explicit CudaArray(std::size_t size);
    ~CudaArray();

    operator T*();

    std::size_t size() const;

    void copy(const T* source);

    template<typename SourceType>
    void upload(const SourceType* source);

    template<typename DestType>
    void download(DestType* destination);

    void clear();

    static std::size_t estimate_mem(std::size_t count);

private:
    void allocate();
    void release();

    T* data = nullptr;
    std::size_t element_count = 0;
};

} // namespace cuda
} // namespace fem

#include "cuda_array.tpp"

#endif
