#include "solve_sparse_indirect.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../cuda/assert_cuda.h"
#include "../../cuda/cuda_array.h"
#include "../../cuda/cuda_csr.h"
#include "../../cuda/cuda_defs.h"
#include "../../cuda/cuda_vec.h"

#include <Eigen/OrderingMethods>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace fem::solver::detail {

#ifdef SUPPORT_GPU
namespace {

constexpr CudaPrecision BENCHMARK_TOLERANCE = static_cast<CudaPrecision>(1e-12);
constexpr Precision VALIDATION_TOLERANCE = static_cast<Precision>(1e-8);
constexpr int BENCHMARK_MAX_ITERATIONS = 100000;
constexpr int LOG_INTERVAL = 1000;

enum class Ordering {
    NATURAL,
    AMD,
    COLAMD
};

struct Variant {
    const char* name;
    Ordering ordering;
    bool scale;
    bool ic0;
    Precision shift;
    cusparseSpMVAlg_t spmv_alg;
};

struct Result {
    std::string name;
    bool converged = false;
    bool valid = false;
    int iterations = 0;
    Precision recursive_residual = std::numeric_limits<Precision>::infinity();
    Precision true_residual = std::numeric_limits<Precision>::infinity();
    double transform_ms = 0.0;
    double setup_ms = 0.0;
    double solve_ms = 0.0;
    double total_ms = 0.0;
    std::string error;
};

struct Prepared {
    using Permutation = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;

    SparseMatrix A;
    DynamicMatrix b;
    DynamicVector D;
    Permutation P;
    bool scaled = false;
    bool permuted = false;

    explicit Prepared(Eigen::Index n)
        : D(DynamicVector::Ones(n))
        , P(n) {
        P.setIdentity();
    }
};

struct SolveResult {
    DynamicMatrix x;
    bool converged = true;
    int iterations = 0;
    Precision residual = 0;
    double setup_ms = 0.0;
    double solve_ms = 0.0;
};

struct DiagonalStats {
    Precision min = std::numeric_limits<Precision>::infinity();
    Precision max = -std::numeric_limits<Precision>::infinity();
    Eigen::Index missing = 0;
    Eigen::Index nonpositive = 0;
    Eigen::Index nonfinite = 0;
};

double elapsed_ms(const std::chrono::steady_clock::time_point& begin) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
}

DiagonalStats diagonal_stats(const SparseMatrix& A) {
    DiagonalStats stats;

    for (Eigen::Index i = 0; i < A.rows(); ++i) {
        const Precision aii = A.coeff(i, i);

        if (aii == Precision{0}) {
            ++stats.missing;
        }
        if (!std::isfinite(aii)) {
            ++stats.nonfinite;
            continue;
        }
        if (aii <= Precision{0}) {
            ++stats.nonpositive;
        }

        stats.min = std::min(stats.min, aii);
        stats.max = std::max(stats.max, aii);
    }

    return stats;
}

void log_diagonal_stats(const SparseMatrix& A, const char* label) {
    const DiagonalStats stats = diagonal_stats(A);
    const Precision ratio =
        stats.min > Precision{0}
            ? stats.max / stats.min
            : std::numeric_limits<Precision>::infinity();

    logging::info(true,
                  "    ", label,
                  " diagonal: min=", stats.min,
                  ", max=", stats.max,
                  ", max/min=", ratio,
                  ", zero=", stats.missing,
                  ", nonpositive=", stats.nonpositive,
                  ", nonfinite=", stats.nonfinite);
}

Eigen::Index original_active_row(Eigen::Index prepared_row,
                                 const Prepared& prepared) {
    if (!prepared.permuted) {
        return prepared_row;
    }
    return prepared.P.indices()(prepared_row);
}

void factor_ic0(cuda::CudaCSR& mat,
                const SparseMatrix& host_matrix,
                const Prepared& prepared) {
    cusparseMatDescr_t descr;
    csric02Info_t info;
    int buffer_size = 0;

    runtime_check_cuda(cusparseCreateCsric02Info(&info));
    runtime_check_cuda(cusparseCreateMatDescr(&descr));
    runtime_check_cuda(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    runtime_check_cuda(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

    runtime_check_cuda(CUSOLV_CSRIC_BUF(cuda::manager.handle_cusparse,
                                        mat.cols(), mat.nnz(), descr,
                                        mat.val_ptr(), mat.row_ptr(), mat.col_ind(),
                                        info, &buffer_size));

    cuda::CudaArray<char> buffer{std::max<size_t>(static_cast<size_t>(buffer_size), 1)};

    runtime_check_cuda(CUSOLV_CSRIC_ANA(cuda::manager.handle_cusparse,
                                        mat.cols(), mat.nnz(), descr,
                                        mat.val_ptr(), mat.row_ptr(), mat.col_ind(),
                                        info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));

    int zero_pivot = -1;
    cusparseStatus_t zero_pivot_status =
        cusparseXcsric02_zeroPivot(cuda::manager.handle_cusparse, info, &zero_pivot);

    if (zero_pivot_status == CUSPARSE_STATUS_ZERO_PIVOT) {
        const Eigen::Index row = static_cast<Eigen::Index>(zero_pivot);
        const Eigen::Index original_row = original_active_row(row, prepared);
        const Precision aii = host_matrix.coeff(row, row);
        logging::error(false,
                       "IC(0) analysis failed: structural zero pivot at prepared row ", row,
                       ", original active row ", original_row,
                       ", diagonal=", aii);
    }
    runtime_check_cuda(zero_pivot_status);

    runtime_check_cuda(CUSOLV_CSRIC(cuda::manager.handle_cusparse,
                                    mat.cols(), mat.nnz(), descr,
                                    mat.val_ptr(), mat.row_ptr(), mat.col_ind(),
                                    info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));

    zero_pivot = -1;
    zero_pivot_status =
        cusparseXcsric02_zeroPivot(cuda::manager.handle_cusparse, info, &zero_pivot);

    if (zero_pivot_status == CUSPARSE_STATUS_ZERO_PIVOT) {
        const Eigen::Index row = static_cast<Eigen::Index>(zero_pivot);
        const Eigen::Index original_row = original_active_row(row, prepared);
        const Precision aii = host_matrix.coeff(row, row);
        logging::error(false,
                       "IC(0) factorization failed: numerical zero pivot at prepared row ", row,
                       ", original active row ", original_row,
                       ", diagonal=", aii);
    }
    runtime_check_cuda(zero_pivot_status);

    runtime_check_cuda(cusparseDestroyCsric02Info(info));
    runtime_check_cuda(cusparseDestroyMatDescr(descr));
}

DynamicVector diagonal_scaling(const SparseMatrix& A) {
    DynamicVector D(A.rows());

    for (Eigen::Index i = 0; i < A.rows(); ++i) {
        const Precision aii = A.coeff(i, i);
        logging::error(std::isfinite(aii) && aii > Precision{0},
                       "Diagonal scaling requires positive finite diagonal at row ", i,
                       ", got ", aii);
        D(i) = Precision{1} / std::sqrt(aii);
    }

    return D;
}

void scale_matrix(SparseMatrix& A, const DynamicVector& D) {
    for (Eigen::Index outer = 0; outer < A.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator it(A, outer); it; ++it) {
            it.valueRef() *= D(it.row()) * D(it.col());
        }
    }
}

void shift_preconditioner(SparseMatrix& A, Precision shift) {
    if (shift == Precision{0}) {
        return;
    }

    for (Eigen::Index i = 0; i < A.rows(); ++i) {
        const Precision aii = A.coeff(i, i);
        logging::error(std::isfinite(aii) && aii > Precision{0},
                       "IC(0) shift requires positive finite diagonal at row ", i,
                       ", got ", aii);
        A.coeffRef(i, i) += shift * aii;
    }

    A.makeCompressed();
}

Prepared prepare(const SparseMatrix& A,
                 const DynamicMatrix& b,
                 const Variant& variant) {
    Prepared out(A.rows());

    if (variant.ordering == Ordering::AMD) {
        Eigen::AMDOrdering<int> ordering;
        ordering(A.selfadjointView<Eigen::Lower>(), out.P);
        out.permuted = true;
    } else if (variant.ordering == Ordering::COLAMD) {
        Eigen::COLAMDOrdering<int> ordering;
        ordering(A, out.P);
        out.permuted = true;
    }

    if (out.permuted) {
        out.A = out.P * A * out.P.transpose();
        out.b = out.P * b;
    } else {
        out.A = A;
        out.b = b;
    }

    if (variant.scale) {
        out.D = diagonal_scaling(out.A);
        scale_matrix(out.A, out.D);

        for (Eigen::Index col = 0; col < out.b.cols(); ++col) {
            out.b.col(col).array() *= out.D.array();
        }

        out.scaled = true;
    }

    out.A.makeCompressed();
    return out;
}

DynamicMatrix restore(DynamicMatrix x, const Prepared& prepared) {
    if (prepared.scaled) {
        for (Eigen::Index col = 0; col < x.cols(); ++col) {
            x.col(col).array() *= prepared.D.array();
        }
    }

    if (prepared.permuted) {
        x = (prepared.P.transpose() * x).eval();
    }

    return x;
}

Precision true_residual(const SparseMatrix& A,
                        const DynamicMatrix& b,
                        const DynamicMatrix& x) {
    Precision worst = 0;

    for (Eigen::Index col = 0; col < b.cols(); ++col) {
        const Precision norm_b = b.col(col).norm();
        if (norm_b == Precision{0}) {
            continue;
        }

        const DynamicVector r = A * x.col(col) - b.col(col);
        worst = std::max(worst, r.norm() / norm_b);
    }

    return worst;
}

SolveResult solve_variant(Prepared& prepared,
                          const Variant& variant) {
    SparseMatrix& A = prepared.A;
    DynamicMatrix& b = prepared.b;

    const auto N = A.cols();
    const auto nnz = A.nonZeros();

    SolveResult result;
    result.x = DynamicMatrix::Zero(N, b.cols());

    const auto setup_begin = std::chrono::steady_clock::now();

    std::unique_ptr<cuda::CudaCSR> prec;
    std::unique_ptr<cuda::CudaCSR> mat_gpu;

    if (variant.ic0) {
        SparseMatrix P = A;
        shift_preconditioner(P, variant.shift);
        log_diagonal_stats(P, variant.name);

        prec = std::make_unique<cuda::CudaCSR>(P);
        mat_gpu = std::make_unique<cuda::CudaCSR>(A, *prec);
        factor_ic0(*prec, P, prepared);
    } else {
        mat_gpu = std::make_unique<cuda::CudaCSR>(A);
    }

    cuda::CudaVector x{static_cast<int>(N)};
    cuda::CudaVector r{static_cast<int>(N)};
    cuda::CudaVector z{static_cast<int>(N)};
    cuda::CudaVector p{static_cast<int>(N)};
    cuda::CudaVector Ap{static_cast<int>(N)};
    cuda::CudaVector tmp{static_cast<int>(N)};

    cusparseSpMatDescr_t descr_A{};
    runtime_check_cuda(cusparseCreateCsr(&descr_A, N, N, nnz,
                                         mat_gpu->row_ptr(),
                                         mat_gpu->col_ind(),
                                         mat_gpu->val_ptr(),
                                         CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_BASE_ZERO,
                                         CUDA_P_TYPE));

    cusparseSpMatDescr_t descr_L{};
    cusparseSpSVDescr_t spsv_fwd{};
    cusparseSpSVDescr_t spsv_bwd{};

    if (variant.ic0) {
        runtime_check_cuda(cusparseCreateCsr(&descr_L, N, N, nnz,
                                             prec->row_ptr(),
                                             prec->col_ind(),
                                             prec->val_ptr(),
                                             CUSPARSE_INDEX_32I,
                                             CUSPARSE_INDEX_32I,
                                             CUSPARSE_INDEX_BASE_ZERO,
                                             CUDA_P_TYPE));

        auto fill_mode = CUSPARSE_FILL_MODE_LOWER;
        runtime_check_cuda(cusparseSpMatSetAttribute(
            descr_L, CUSPARSE_SPMAT_FILL_MODE, &fill_mode, sizeof(fill_mode)));

        auto diag_type = CUSPARSE_DIAG_TYPE_NON_UNIT;
        runtime_check_cuda(cusparseSpMatSetAttribute(
            descr_L, CUSPARSE_SPMAT_DIAG_TYPE, &diag_type, sizeof(diag_type)));

        runtime_check_cuda(cusparseSpSV_createDescr(&spsv_fwd));
        runtime_check_cuda(cusparseSpSV_createDescr(&spsv_bwd));
    }

    CudaPrecision one = 1;
    CudaPrecision zero = 0;

    size_t spmv_buffer_size = 0;
    runtime_check_cuda(cusparseSpMV_bufferSize(
        cuda::manager.handle_cusparse,
        CUSPARSE_OPERATION_NON_TRANSPOSE,
        &one, descr_A, p, &zero, Ap,
        CUDA_P_TYPE, variant.spmv_alg, &spmv_buffer_size));

    cuda::CudaArray<char> spmv_buffer{std::max<size_t>(spmv_buffer_size, 1)};

    std::unique_ptr<cuda::CudaArray<char>> fwd_buffer;
    std::unique_ptr<cuda::CudaArray<char>> bwd_buffer;

    if (variant.ic0) {
        size_t fwd_size = 0;
        size_t bwd_size = 0;

        runtime_check_cuda(cusparseSpSV_bufferSize(
            cuda::manager.handle_cusparse,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one, descr_L, r, tmp,
            CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_fwd, &fwd_size));

        runtime_check_cuda(cusparseSpSV_bufferSize(
            cuda::manager.handle_cusparse,
            CUSPARSE_OPERATION_TRANSPOSE,
            &one, descr_L, tmp, z,
            CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_bwd, &bwd_size));

        fwd_buffer = std::make_unique<cuda::CudaArray<char>>(std::max<size_t>(fwd_size, 1));
        bwd_buffer = std::make_unique<cuda::CudaArray<char>>(std::max<size_t>(bwd_size, 1));

        runtime_check_cuda(cusparseSpSV_analysis(
            cuda::manager.handle_cusparse,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one, descr_L, r, tmp,
            CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_fwd, *fwd_buffer));

        runtime_check_cuda(cusparseSpSV_analysis(
            cuda::manager.handle_cusparse,
            CUSPARSE_OPERATION_TRANSPOSE,
            &one, descr_L, tmp, z,
            CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_bwd, *bwd_buffer));
    }

    runtime_check_cuda(cudaDeviceSynchronize());
    result.setup_ms = elapsed_ms(setup_begin);

    const auto solve_begin = std::chrono::steady_clock::now();
    const int max_iterations =
        std::min(static_cast<int>(N), BENCHMARK_MAX_ITERATIONS);

    for (Eigen::Index col = 0; col < b.cols(); ++col) {
        if (b.col(col).isZero()) {
            continue;
        }

        x.clear();
        r.upload(b.col(col).data());

        if (variant.ic0) {
            runtime_check_cuda(cusparseSpSV_solve(
                cuda::manager.handle_cusparse,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                &one, descr_L, r, tmp,
                CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_fwd));

            runtime_check_cuda(cusparseSpSV_solve(
                cuda::manager.handle_cusparse,
                CUSPARSE_OPERATION_TRANSPOSE,
                &one, descr_L, tmp, z,
                CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_bwd));

            p.copy(z);
        } else {
            p.copy(r);
        }

        const CudaPrecision norm_b =
            static_cast<CudaPrecision>(b.col(col).norm());

        CudaPrecision rel_res = 1;
        int k = 0;

        for (k = 1; k <= max_iterations; ++k) {
            CudaPrecision rz = 0;
            CudaPrecision pAp = 0;

            if (variant.ic0) {
                runtime_check_cuda(CUBLAS_DOT(
                    cuda::manager.handle_cublas, N, r, 1, z, 1, &rz));
            } else {
                runtime_check_cuda(CUBLAS_DOT(
                    cuda::manager.handle_cublas, N, r, 1, r, 1, &rz));
            }

            runtime_check_cuda(cusparseSpMV(
                cuda::manager.handle_cusparse,
                CUSPARSE_OPERATION_NON_TRANSPOSE,
                &one, descr_A, p, &zero, Ap,
                CUDA_P_TYPE, variant.spmv_alg, spmv_buffer));

            runtime_check_cuda(CUBLAS_DOT(
                cuda::manager.handle_cublas, N, Ap, 1, p, 1, &pAp));

            if (!std::isfinite(rz) || !std::isfinite(pAp) || pAp <= CudaPrecision{0}) {
                rel_res = std::numeric_limits<CudaPrecision>::infinity();
                break;
            }

            CudaPrecision alpha = rz / pAp;
            CudaPrecision neg_alpha = -alpha;

            runtime_check_cuda(CUBLAS_AXPY(
                cuda::manager.handle_cublas, N, &alpha, p, 1, x, 1));

            runtime_check_cuda(CUBLAS_AXPY(
                cuda::manager.handle_cublas, N, &neg_alpha, Ap, 1, r, 1));

            CudaPrecision r_norm = 0;
            runtime_check_cuda(CUBLAS_NRM(
                cuda::manager.handle_cublas, N, r, 1, &r_norm));

            rel_res = r_norm / norm_b;

            if (k % LOG_INTERVAL == 0) {
                logging::info(true,
                              "    ", variant.name,
                              " RHS ", col,
                              " iteration ", k,
                              " relative residual: ", rel_res);
            }

            if (!std::isfinite(rel_res) || rel_res < BENCHMARK_TOLERANCE) {
                break;
            }

            CudaPrecision rz_new = 0;

            if (variant.ic0) {
                runtime_check_cuda(cusparseSpSV_solve(
                    cuda::manager.handle_cusparse,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    &one, descr_L, r, tmp,
                    CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_fwd));

                runtime_check_cuda(cusparseSpSV_solve(
                    cuda::manager.handle_cusparse,
                    CUSPARSE_OPERATION_TRANSPOSE,
                    &one, descr_L, tmp, z,
                    CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_bwd));

                runtime_check_cuda(CUBLAS_DOT(
                    cuda::manager.handle_cublas, N, r, 1, z, 1, &rz_new));
            } else {
                runtime_check_cuda(CUBLAS_DOT(
                    cuda::manager.handle_cublas, N, r, 1, r, 1, &rz_new));
            }

            if (!std::isfinite(rz_new) || rz == CudaPrecision{0}) {
                rel_res = std::numeric_limits<CudaPrecision>::infinity();
                break;
            }

            CudaPrecision beta = rz_new / rz;

            runtime_check_cuda(CUBLAS_SCAL(
                cuda::manager.handle_cublas, N, &beta, p, 1));

            cuda::CudaVector& direction = variant.ic0 ? z : r;
            runtime_check_cuda(CUBLAS_AXPY(
                cuda::manager.handle_cublas, N, &one, direction, 1, p, 1));
        }

        x.download(result.x.col(col).data());

        result.iterations =
            std::max(result.iterations, std::min(k, max_iterations));
        result.residual =
            std::max(result.residual, static_cast<Precision>(rel_res));

        if (!std::isfinite(rel_res) || rel_res >= BENCHMARK_TOLERANCE) {
            result.converged = false;
        }
    }

    runtime_check_cuda(cudaDeviceSynchronize());
    result.solve_ms = elapsed_ms(solve_begin);

    if (variant.ic0) {
        runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_fwd));
        runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_bwd));
        runtime_check_cuda(cusparseDestroySpMat(descr_L));
    }

    runtime_check_cuda(cusparseDestroySpMat(descr_A));

    return result;
}

void print_summary(std::vector<Result> results) {
    std::sort(results.begin(), results.end(),
              [](const Result& a, const Result& b) {
                  if (a.valid != b.valid) {
                      return a.valid > b.valid;
                  }
                  return a.total_ms < b.total_ms;
              });

    logging::info(true, "");
    logging::info(true, "============================================================================================================");
    logging::info(true, "GPU INDIRECT SOLVER BENCHMARK SUMMARY");
    logging::info(true, "============================================================================================================");
    logging::info(true,
                  std::setw(36), std::left, "method",
                  std::setw(8), std::right, "valid",
                  std::setw(10), std::right, "iters",
                  std::setw(14), std::right, "transform",
                  std::setw(14), std::right, "setup",
                  std::setw(14), std::right, "solve",
                  std::setw(14), std::right, "total",
                  std::setw(18), std::right, "true residual");

    for (const Result& r : results) {
        logging::info(true,
                      std::setw(36), std::left, r.name,
                      std::setw(8), std::right, (r.valid ? "yes" : "no"),
                      std::setw(10), std::right, r.iterations,
                      std::setw(14), std::right, r.transform_ms,
                      std::setw(14), std::right, r.setup_ms,
                      std::setw(14), std::right, r.solve_ms,
                      std::setw(14), std::right, r.total_ms,
                      std::setw(18), std::right, r.true_residual);

        if (!r.error.empty()) {
            logging::info(true, "    error: ", r.error);
        }
    }

    logging::info(true, "============================================================================================================");
}

} // namespace
#endif

DynamicMatrix solve_indirect_gpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs) {
#ifndef SUPPORT_GPU
    logging::info(true,
                  "This build does not support gpu-accelerated solving, falling back to cpu");
    return solve_indirect_cpu(mat, rhs);
#else
    logging::error(mat.rows() == mat.cols(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "rhs row count must match matrix size");

    cuda::manager.create_cuda();

    const std::vector<Variant> variants = {
        {"01 CG natural / ALG1",                 Ordering::NATURAL, false, false, Precision{0}, CUSPARSE_SPMV_CSR_ALG1},
        {"02 diagonal-scaled CG",                Ordering::NATURAL, true,  false, Precision{0}, CUSPARSE_SPMV_CSR_ALG1},
        {"03 AMD + IC0",                         Ordering::AMD,     false, true,  Precision{0}, CUSPARSE_SPMV_CSR_ALG1},
        {"04 AMD + scaled + IC0",                Ordering::AMD,     true,  true,  Precision{0}, CUSPARSE_SPMV_CSR_ALG1},
        {"05 AMD + scaled + IC0 / shift 1e-2",  Ordering::AMD,     true,  true,  static_cast<Precision>(1e-2), CUSPARSE_SPMV_CSR_ALG1},
        {"06 scaled + IC0",                      Ordering::NATURAL, true,  true,  Precision{0}, CUSPARSE_SPMV_CSR_ALG1},
        {"07 scaled + IC0 / shift 1e-4",        Ordering::NATURAL, true,  true,  static_cast<Precision>(1e-4), CUSPARSE_SPMV_CSR_ALG1},
        {"08 scaled + IC0 / shift 1e-2",        Ordering::NATURAL, true,  true,  static_cast<Precision>(1e-2), CUSPARSE_SPMV_CSR_ALG1},
        {"09 scaled + IC0 / shift 1e-1",        Ordering::NATURAL, true,  true,  static_cast<Precision>(1e-1), CUSPARSE_SPMV_CSR_ALG1},
        {"10 scaled + IC0 / shift 1",           Ordering::NATURAL, true,  true,  Precision{1}, CUSPARSE_SPMV_CSR_ALG1},
        {"11 IC0 natural / shift 1e-1",         Ordering::NATURAL, false, true,  static_cast<Precision>(1e-1), CUSPARSE_SPMV_CSR_ALG1},
        {"12 IC0 natural / shift 1",            Ordering::NATURAL, false, true,  Precision{1}, CUSPARSE_SPMV_CSR_ALG1},
    };

    logging::info(true, "");
    logging::info(true, "============================================================================================================");
    logging::info(true, "GPU INDIRECT SOLVER BENCHMARK");
    logging::info(true, "============================================================================================================");
    logging::info(true, "N                    : ", mat.rows());
    logging::info(true, "nnz                  : ", mat.nonZeros());
    logging::info(true, "RHS columns          : ", rhs.cols());
    logging::info(true, "PCG tolerance        : ", BENCHMARK_TOLERANCE);
    logging::info(true, "validation tolerance : ", VALIDATION_TOLERANCE);
    logging::info(true, "max iterations       : ", BENCHMARK_MAX_ITERATIONS);
    logging::info(true, "variants             : ", variants.size());
    log_diagonal_stats(mat, "original matrix");
    logging::info(true, "============================================================================================================");

    std::vector<Result> results;
    DynamicMatrix best_solution;
    DynamicMatrix fallback_solution;
    bool have_valid = false;
    double best_time = std::numeric_limits<double>::infinity();
    Precision best_residual = std::numeric_limits<Precision>::infinity();

    for (size_t i = 0; i < variants.size(); ++i) {
        const Variant& variant = variants[i];
        Result result;
        result.name = variant.name;

        logging::info(true, "");
        logging::info(true, "[", i + 1, "/", variants.size(), "] ", variant.name);

        const auto total_begin = std::chrono::steady_clock::now();

        try {
            const auto transform_begin = std::chrono::steady_clock::now();
            Prepared prepared = prepare(mat, rhs, variant);
            result.transform_ms = elapsed_ms(transform_begin);

            SolveResult solve = solve_variant(prepared, variant);
            result.setup_ms = solve.setup_ms;
            result.solve_ms = solve.solve_ms;
            result.iterations = solve.iterations;
            result.recursive_residual = solve.residual;
            result.converged = solve.converged;

            DynamicMatrix x = restore(std::move(solve.x), prepared);
            result.true_residual = true_residual(mat, rhs, x);
            result.valid = result.converged &&
                           std::isfinite(result.true_residual) &&
                           result.true_residual <= VALIDATION_TOLERANCE;
            result.total_ms = elapsed_ms(total_begin);

            logging::info(true,
                          "    iterations=", result.iterations,
                          ", recursive residual=", result.recursive_residual,
                          ", true residual=", result.true_residual,
                          ", transform=", result.transform_ms, " ms",
                          ", setup=", result.setup_ms, " ms",
                          ", solve=", result.solve_ms, " ms",
                          ", total=", result.total_ms, " ms",
                          ", valid=", (result.valid ? "yes" : "no"));

            if (std::isfinite(result.true_residual) &&
                result.true_residual < best_residual) {
                best_residual = result.true_residual;
                fallback_solution = x;
            }

            if (result.valid && result.total_ms < best_time) {
                best_time = result.total_ms;
                best_solution = std::move(x);
                have_valid = true;
            }
        } catch (const std::exception& e) {
            result.error = e.what();
            result.total_ms = elapsed_ms(total_begin);
            logging::info(true, "    FAILED: ", result.error);
        } catch (...) {
            result.error = "unknown exception";
            result.total_ms = elapsed_ms(total_begin);
            logging::info(true, "    FAILED: unknown exception");
        }

        results.push_back(std::move(result));

        cudaGetLastError();
        cudaDeviceSynchronize();
        cudaGetLastError();
    }

    print_summary(results);

    if (have_valid) {
        logging::info(true, "Returning fastest valid benchmark solution");
        return best_solution;
    }

    logging::warning(true,
                     "No variant reached validation tolerance; "
                     "returning solution with smallest true residual");
    return fallback_solution;
#endif
}

} // namespace fem::solver::detail