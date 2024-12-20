#pragma once

#ifdef SUPPORT_GPU

#include "../core/types_eig.h"
#include "cuda_array.h"

#include <memory>

namespace cuda{

struct CudaCSR{

    private:
    // create csr arrays on the gpu
    std::shared_ptr<CudaArray<CudaPrecision>> m_val_ptr;
    std::shared_ptr<CudaArray<int          >> m_col_ind;
    std::shared_ptr<CudaArray<int          >> m_row_ptr;

    // other data
    size_t m_nnz;
    size_t m_cols;

    public:
    CudaCSR(SparseMatrix &matrix)
        : m_val_ptr(new CudaArray<CudaPrecision>(matrix.nonZeros()))
        , m_col_ind(new CudaArray<int          >(matrix.nonZeros()))
        , m_row_ptr(new CudaArray<int          >(matrix.rows() + 1))
        , m_nnz(matrix.nonZeros())
        , m_cols(matrix.cols()){
        m_val_ptr->upload(matrix.valuePtr());
        m_col_ind->upload(matrix.innerIndexPtr());
        m_row_ptr->upload(matrix.outerIndexPtr());
    }

    CudaCSR(SparseMatrix &matrix, CudaCSR &similar)
            : m_val_ptr(new CudaArray<CudaPrecision>(matrix.nonZeros()))
            , m_col_ind(similar.m_col_ind)
            , m_row_ptr(similar.m_row_ptr)
            , m_nnz(similar.m_nnz)
            , m_cols(similar.m_cols){
        runtime_check(matrix.nonZeros() == m_nnz, "cannot construct matrix with same column indices and row pointers");
        runtime_check(matrix.cols() == m_cols, "cannot construct matrix with same column indices and row pointers");
        m_val_ptr->upload(matrix.valuePtr());
    }

    void download(SparseMatrix &matrix){
        runtime_check(matrix.nonZeros() == m_nnz, "cannot construct matrix with same column indices and row pointers");
        runtime_check(matrix.cols() == m_cols, "cannot construct matrix with same column indices and row pointers");

        m_val_ptr->download(matrix.valuePtr());
        m_col_ind->download(matrix.innerIndexPtr());
        m_row_ptr->download(matrix.outerIndexPtr());
    }

    CudaArray<CudaPrecision>& val_ptr(){
        return *m_val_ptr;
    }
    CudaArray<int      >& col_ind(){
        return *m_col_ind;
    }
    CudaArray<int      >& row_ptr(){
        return *m_row_ptr;
    }
    size_t nnz() const {
        return m_nnz;
    }
    size_t cols() const {
        return m_cols;
    }

    static size_t estimate_mem(SparseMatrix &matrix, bool only_values=false){
        size_t s = CudaArray<CudaPrecision>::estimate_mem(matrix.nonZeros());
        if(!only_values){
            s += CudaArray<int>::estimate_mem(matrix.nonZeros());
            s += CudaArray<int>::estimate_mem(matrix.rows() + 1);
        }
        return s;
    }

};

}

#endif