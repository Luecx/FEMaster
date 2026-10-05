/**
 * @file cuda_csr.h
 * @brief Defines device CSR storage for sparse GPU operations.
 *
 * CudaCSR uploads the CSC representation of the transposed Eigen matrix,
 * which has the CSR layout required for the original matrix on the GPU.
 * Values use CudaPrecision while row offsets and column indices use int.
 * CudaArray owns device allocations; solver descriptors and factorization
 * remain responsibilities of the calling GPU solver.
 *
 * @see CudaCSR
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#ifdef SUPPORT_GPU

#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "cuda_array.h"

#include <memory>

namespace fem::cuda{
/**
 * Owns device values and shared CSR topology for an Eigen sparse matrix.
 *
 * Upload constructors allocate values independently. Construction from a similar
 * instance shares its row offsets and column indices through shared ownership;
 * ordinary copies share all device arrays.
 * Reusing topology requires the same sparsity pattern, not merely matching sizes.
 * Upload and download transpose the column-major Eigen representation to obtain
 * row-major CSR storage. Downloads may replace the destination topology while
 * requiring matching nonzero and column counts. Solver state is kept externally.
 */
struct CudaCSR{
    private:
    // Device values and shared 32-bit CSR topology.
    std::shared_ptr<CudaArray<CudaPrecision>> m_val_ptr;
    std::shared_ptr<CudaArray<int          >> m_col_ind;
    std::shared_ptr<CudaArray<int          >> m_row_ptr;

    // Nonzero and column counts used to validate compatible transfers.
    size_t m_nnz;
    size_t m_cols;

    public:
    // Upload independent storage, or reuse the topology of a matching matrix.
    CudaCSR(SparseMatrix &matrix)
        : m_val_ptr(new CudaArray<CudaPrecision>(matrix.nonZeros()))
        , m_col_ind(new CudaArray<int          >(matrix.nonZeros()))
        , m_row_ptr(new CudaArray<int          >(matrix.rows() + 1))
        , m_nnz(matrix.nonZeros())
        , m_cols(matrix.cols()) {
        // CSC storage of A^T has the same layout as CSR storage of A
        SparseMatrix matrix_t = matrix.transpose();
        matrix_t.makeCompressed();

        m_val_ptr->upload(matrix_t.valuePtr());
        m_col_ind->upload(matrix_t.innerIndexPtr());
        m_row_ptr->upload(matrix_t.outerIndexPtr());
    }

    CudaCSR(SparseMatrix &matrix, CudaCSR &similar)
            : m_val_ptr(new CudaArray<CudaPrecision>(matrix.nonZeros()))
            , m_col_ind(similar.m_col_ind)
            , m_row_ptr(similar.m_row_ptr)
            , m_nnz(similar.m_nnz)
            , m_cols(similar.m_cols) {
        // Eigen dimensions are nonnegative; compare in the device count type.
        runtime_check(static_cast<size_t>(matrix.nonZeros()) == m_nnz,
            "cannot construct matrix with same column indices and row pointers");
        runtime_check(static_cast<size_t>(matrix.cols()) == m_cols,
            "cannot construct matrix with same column indices and row pointers");

        // Match the shared CSR ordering through compressed CSC storage of A^T.
        SparseMatrix matrix_t = matrix.transpose();
        matrix_t.makeCompressed();
        m_val_ptr->upload(matrix_t.valuePtr());
    }

    // Download into a matrix with compatible nonzero and column counts.
    void download(SparseMatrix &matrix) {
        // Eigen dimensions are nonnegative; compare in the device count type.
        runtime_check(static_cast<size_t>(matrix.nonZeros()) == m_nnz,
            "cannot construct matrix with same column indices and row pointers");
        runtime_check(static_cast<size_t>(matrix.cols()) == m_cols,
            "cannot construct matrix with same column indices and row pointers");

        // Restore device CSR arrays through CSC storage of the transpose.
        SparseMatrix matrix_t = matrix.transpose();
        matrix_t.makeCompressed();

        m_val_ptr->download(matrix_t.valuePtr());
        m_col_ind->download(matrix_t.innerIndexPtr());
        m_row_ptr->download(matrix_t.outerIndexPtr());
        matrix = matrix_t.transpose();
    }

    CudaArray<CudaPrecision>& val_ptr() {
        return *m_val_ptr;
    }
    CudaArray<int      >& col_ind() {
        return *m_col_ind;
    }
    CudaArray<int      >& row_ptr() {
        return *m_row_ptr;
    }
    size_t nnz() const {
        return m_nnz;
    }
    size_t cols() const {
        return m_cols;
    }

    static size_t estimate_mem(SparseMatrix &matrix, bool only_values=false) {
        size_t s = CudaArray<CudaPrecision>::estimate_mem(matrix.nonZeros());
        if(!only_values) {
            s += CudaArray<int>::estimate_mem(matrix.nonZeros());
            s += CudaArray<int>::estimate_mem(matrix.rows() + 1);
        }
        return s;
    }
};
}

#endif
