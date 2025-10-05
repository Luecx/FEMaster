/**
 * @file cuda_array.tpp
 * @brief Template implementation for `CudaArray`.
 */

#ifdef SUPPORT_GPU

namespace fem {
namespace cuda {

template<typename T>
CudaArray<T>::CudaArray(std::size_t size) : element_count(size) {
    runtime_check(size > 0, "size cannot be less than 1");
    allocate();
}

template<typename T>
CudaArray<T>::~CudaArray() {
    release();
}

template<typename T>
CudaArray<T>::operator T*() {
    return data;
}

template<typename T>
std::size_t CudaArray<T>::size() const {
    return element_count;
}

template<typename T>
void CudaArray<T>::copy(const T* source) {
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    runtime_check_cuda(cudaMemcpy(data, source, size() * sizeof(T), cudaMemcpyDeviceToDevice));
}

template<typename T>
template<typename SourceType>
void CudaArray<T>::upload(const SourceType* source) {
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    if constexpr (std::is_same_v<SourceType, T>) {
        runtime_check_cuda(cudaMemcpy(data, source, size() * sizeof(T), cudaMemcpyHostToDevice));
    } else {
        std::vector<T> intermediate(size());
        for (std::size_t i = 0; i < size(); ++i) {
            intermediate[i] = static_cast<T>(source[i]);
        }
        runtime_check_cuda(cudaMemcpy(data, intermediate.data(), size() * sizeof(T), cudaMemcpyHostToDevice));
    }
}

template<typename T>
template<typename DestType>
void CudaArray<T>::download(DestType* destination) {
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    if constexpr (std::is_same_v<DestType, T>) {
        runtime_check_cuda(cudaMemcpy(destination, data, size() * sizeof(T), cudaMemcpyDeviceToHost));
    } else {
        std::vector<T> intermediate(size());
        runtime_check_cuda(cudaMemcpy(intermediate.data(), data, size() * sizeof(T), cudaMemcpyDeviceToHost));
        for (std::size_t i = 0; i < size(); ++i) {
            destination[i] = static_cast<DestType>(intermediate[i]);
        }
    }
}

template<typename T>
void CudaArray<T>::clear() {
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    runtime_check_cuda(cudaMemset(data, 0, size() * sizeof(T)));
}

template<typename T>
std::size_t CudaArray<T>::estimate_mem(std::size_t count) {
    return sizeof(T) * count;
}

template<typename T>
void CudaArray<T>::allocate() {
    if (data != nullptr) {
        return;
    }
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    runtime_check_cuda(cudaMalloc(reinterpret_cast<void**>(&data), size() * sizeof(T)));
    clear();
}

template<typename T>
void CudaArray<T>::release() {
    if (data == nullptr) {
        return;
    }
    runtime_assert(fem::cuda::manager.is_loaded(), "CUDA is not loaded");
    runtime_check_cuda(cudaFree(data));
    data = nullptr;
    element_count = 0;
}

} // namespace cuda
} // namespace fem

#endif
