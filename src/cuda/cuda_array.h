#pragma once

#ifdef SUPPORT_GPU


#include "cuda.h"

namespace cuda {

template<typename T>
class CudaArray {

    protected:
    T* data = nullptr;
    size_t m_size = 0;

    public:
    CudaArray(size_t size) : m_size(size) {
        runtime_check(size > 0, "size cannot be less than 1");
        malloc();
    }

    virtual ~CudaArray() {
        this->free();
    }

    operator T*(){
        return data;
    }

    int size(){
        return m_size;
    }

    private:
    void malloc() {
        if(this->data != nullptr) return;
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaMalloc((void**) &this->data, this->m_size * sizeof(T)));
        clear();
    }

    void free() {
        if(this->data == nullptr) return;
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaFree(this->data));
        this->data = nullptr;
        this->m_size = 0;
    }

    public:
    void copy(const T* source) {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaMemcpy(this->data, source, this->size() * sizeof(T), cudaMemcpyDeviceToDevice));
    }

    template <typename SourceType>
    void upload(const SourceType* source) {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");

        if constexpr (std::is_same_v<SourceType, T>) {
            // If the source is of type T, directly copy
            runtime_check_cuda(cudaMemcpy(this->data, source, this->size() * sizeof(T), cudaMemcpyHostToDevice));
        } else {
            // If the source is not of type T, create an intermediate array of type T
            T* intermediate = new T[this->size()];
            for (size_t i = 0; i < this->size(); i++) {
                intermediate[i] = static_cast<T>(source[i]);
            }

            // Then copy the intermediate array to device
            runtime_check_cuda(cudaMemcpy(this->data, intermediate, this->size() * sizeof(T), cudaMemcpyHostToDevice));

            // Don't forget to delete the intermediate array
            delete[] intermediate;
        }
    }


    template <typename DestType>
    void download(DestType* destination) {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");

        if constexpr (std::is_same_v<DestType, T>) {
            // If the destination is of type T, directly copy
            runtime_check_cuda(cudaMemcpy(destination, this->data, this->size() * sizeof(T), cudaMemcpyDeviceToHost));
        } else {
            // If the destination is not of type T, create an intermediate array of type T
            T* intermediate = new T[this->size()];

            // Copy from device to the intermediate array
            runtime_check_cuda(cudaMemcpy(intermediate, this->data, this->size() * sizeof(T), cudaMemcpyDeviceToHost));

            // Then copy from the intermediate array to the destination
            for (size_t i = 0; i < this->size(); i++) {
                destination[i] = static_cast<DestType>(intermediate[i]);
            }

            // Don't forget to delete the intermediate array
            delete[] intermediate;
        }
    }

    void clear() {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaMemset(this->data, 0, this->size() * sizeof(T)));
    }

    static size_t estimate_mem(size_t s){
        return sizeof(T) * s;
    }
};

};    // namespace array

#endif