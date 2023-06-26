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

    void upload(const T* source) {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaMemcpy(this->data, source, this->size() * sizeof(T), cudaMemcpyHostToDevice));
    }

    void download(T* destination) {
        runtime_assert(cuda::manager.is_loaded(), "CUDA is not loaded");
        runtime_check_cuda(cudaMemcpy(destination, this->data, this->size() * sizeof(T), cudaMemcpyDeviceToHost));
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