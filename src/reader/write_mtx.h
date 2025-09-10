// debug_io.h
#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <type_traits>
#include <Eigen/Sparse>

/**
 * Write an Eigen::SparseMatrix to a simple .mtx file:
 * First line: "rows cols nnz"
 * Then: one line per nonzero: "row col value"
 *
 * @param path      Output file path.
 * @param A         Sparse matrix (need not be compressed).
 * @param one_based If true, write 1-based indices; else 0-based.
 * @param precision Floating output precision (default 17).
 */
template <class Scalar>
inline void write_mtx(const std::string& path,
                      const Eigen::SparseMatrix<Scalar>& A,
                      bool one_based = false,
                      int precision = 17)
{
    std::ofstream os(path, std::ios::out | std::ios::trunc);
    if (!os) throw std::runtime_error("write_mtx: could not open '" + path + "'");

    // We can iterate safely without forcing compressed; but make a local copy compressed if you prefer:
    Eigen::SparseMatrix<Scalar> M = A;
    M.makeCompressed();

    const auto R = static_cast<long long>(M.rows());
    const auto C = static_cast<long long>(M.cols());
    const auto NZ = static_cast<long long>(M.nonZeros());
    os << R << " " << C << " " << NZ << "\n";

    os.setf(std::ios::scientific);
    os << std::setprecision(precision);

    for (int k = 0; k < M.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(M, k); it; ++it)
        {
            long long r = it.row();
            long long c = it.col();
            if (one_based) { ++r; ++c; }
            os << r << " " << c << " " << static_cast<long double>(it.value()) << "\n";
        }
}

/**
 * Convenience overload for dense matrices: internally converts to sparse.
 */
template <class Derived>
inline void write_mtx_dense(const std::string& path,
                            const Eigen::MatrixBase<Derived>& A,
                            bool one_based = false,
                            int precision = 17)
{
    using Scalar = typename Derived::Scalar;
    Eigen::SparseMatrix<Scalar> S = A.sparseView();
    write_mtx(path, S, one_based, precision);
}
