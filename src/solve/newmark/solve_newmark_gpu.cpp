#include "solve_newmark_gpu.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../cuda/assert_cuda.h"
#include "../../cuda/cuda_array.h"
#include "../../cuda/cuda_csr.h"
#include "../../cuda/cuda_defs.h"
#include "../../cuda/cuda_vec.h"

#if defined(SUPPORT_GPU) && defined(USE_CUDSS)
#include <cudss.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace fem::solver::detail {

#ifndef SUPPORT_GPU

NewmarkResult
newmark_linear_gpu(const SparseMatrix&,
                   const SparseMatrix&,
                   const SparseMatrix&,
                   const NewmarkIC&,
                   const NewmarkOpts&,
                   const NewmarkForceBasis&)
{
    logging::error(false, "Newmark GPU path requires SUPPORT_GPU.");
    return {};
}

#elif !defined(USE_CUDSS)

NewmarkResult
newmark_linear_gpu(const SparseMatrix&,
                   const SparseMatrix&,
                   const SparseMatrix&,
                   const NewmarkIC&,
                   const NewmarkOpts&,
                   const NewmarkForceBasis&)
{
    logging::error(false, "Newmark GPU path requires USE_CUDSS.");
    return {};
}

#else

namespace {

constexpr CudaPrecision kOne  = CudaPrecision(1);
constexpr CudaPrecision kZero = CudaPrecision(0);

int checked_size(Eigen::Index n) {
    logging::error(n > 0, "Newmark GPU: system size must be positive.");
    logging::error(n <= std::numeric_limits<int>::max(), "Newmark GPU: system size exceeds cuBLAS int range.");
    return static_cast<int>(n);
}

SparseMatrix compressed_copy(const SparseMatrix& matrix) {
    SparseMatrix copy = matrix;
    copy.makeCompressed();
    return copy;
}

void copy_device_vector(cuda::CudaArray<CudaPrecision>& history,
                        std::size_t offset,
                        cuda::CudaVector& vector,
                        std::size_t n)
{
    runtime_check_cuda(cudaMemcpy(static_cast<CudaPrecision*>(history) + offset,
                                  static_cast<CudaPrecision*>(vector),
                                  n * sizeof(CudaPrecision),
                                  cudaMemcpyDeviceToDevice));
}

const char* cudss_status_name(cudssStatus_t status) {
    switch (status) {
        case CUDSS_STATUS_SUCCESS: return "CUDSS_STATUS_SUCCESS";
        case CUDSS_STATUS_NOT_INITIALIZED: return "CUDSS_STATUS_NOT_INITIALIZED";
        case CUDSS_STATUS_ALLOC_FAILED: return "CUDSS_STATUS_ALLOC_FAILED";
        case CUDSS_STATUS_INVALID_VALUE: return "CUDSS_STATUS_INVALID_VALUE";
        case CUDSS_STATUS_NOT_SUPPORTED: return "CUDSS_STATUS_NOT_SUPPORTED";
        case CUDSS_STATUS_EXECUTION_FAILED: return "CUDSS_STATUS_EXECUTION_FAILED";
        case CUDSS_STATUS_INTERNAL_ERROR: return "CUDSS_STATUS_INTERNAL_ERROR";
        default: return "CUDSS_STATUS_UNKNOWN";
    }
}

void runtime_check_cudss(cudssStatus_t status, const char* call) {
    logging::error(status == CUDSS_STATUS_SUCCESS,
                   "cuDSS call failed: ", call, " status=", cudss_status_name(status),
                   " (", static_cast<int>(status), ")");
}

struct GpuSpmv {
    std::unique_ptr<cuda::CudaCSR> mat;
    cusparseSpMatDescr_t descr{};
    std::unique_ptr<cuda::CudaArray<char>> buffer;
    bool empty = false;

    GpuSpmv(SparseMatrix& matrix, cuda::CudaVector& x, cuda::CudaVector& y)
        : empty(matrix.nonZeros() == 0)
    {
        if (empty) {
            return;
        }

        mat = std::make_unique<cuda::CudaCSR>(matrix);
        runtime_check_cuda(cusparseCreateCsr(&descr,
                                             static_cast<int64_t>(matrix.rows()),
                                             static_cast<int64_t>(matrix.cols()),
                                             static_cast<int64_t>(matrix.nonZeros()),
                                             mat->row_ptr(),
                                             mat->col_ind(),
                                             mat->val_ptr(),
                                             CUSPARSE_INDEX_32I,
                                             CUSPARSE_INDEX_32I,
                                             CUSPARSE_INDEX_BASE_ZERO,
                                             CUDA_P_TYPE));

        std::size_t buffer_size = 0;
        runtime_check_cuda(cusparseSpMV_bufferSize(cuda::manager.handle_cusparse,
                                                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &kOne,
                                                   descr,
                                                   x,
                                                   &kZero,
                                                   y,
                                                   CUDA_P_TYPE,
                                                   CUSPARSE_SPMV_CSR_ALG1,
                                                   &buffer_size));
        buffer = std::make_unique<cuda::CudaArray<char>>(std::max<std::size_t>(1, buffer_size));
    }

    ~GpuSpmv() {
        if (descr) {
            runtime_check_cuda(cusparseDestroySpMat(descr));
        }
    }

    void apply(cuda::CudaVector& x, cuda::CudaVector& y, CudaPrecision alpha, CudaPrecision beta) {
        if (empty) {
            if (beta == kZero) {
                y.clear();
            } else if (beta != kOne) {
                runtime_check_cuda(CUBLAS_SCAL(cuda::manager.handle_cublas,
                                               checked_size(static_cast<Eigen::Index>(y.size())),
                                               &beta,
                                               y,
                                               1));
            }
            return;
        }

        runtime_check_cuda(cusparseSpMV(cuda::manager.handle_cusparse,
                                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                                        &alpha,
                                        descr,
                                        x,
                                        &beta,
                                        y,
                                        CUDA_P_TYPE,
                                        CUSPARSE_SPMV_CSR_ALG1,
                                        *buffer));
    }
};

struct CudssFactorizedSolver {
    int n = 0;
    cuda::CudaCSR mat;

    cudssHandle_t handle {};
    cudssConfig_t config {};
    cudssData_t data {};
    cudssMatrix_t cudss_mat {};
    double factorization_time_ms = 0.0;

    CudssFactorizedSolver(SparseMatrix& matrix, const char* name)
        : n(checked_size(matrix.rows())),
          mat(matrix)
    {
        (void)name;
        logging::error(matrix.rows() == matrix.cols(), "Newmark GPU: cuDSS matrix must be square.");

        Timer timer{};
        timer.start();

        runtime_check_cudss(cudssCreate(&handle), "cudssCreate");
        runtime_check_cudss(cudssConfigCreate(&config), "cudssConfigCreate");
        runtime_check_cudss(cudssDataCreate(handle, &data), "cudssDataCreate");

        runtime_check_cudss(cudssMatrixCreateCsr(&cudss_mat,
                                                 static_cast<int64_t>(matrix.rows()),
                                                 static_cast<int64_t>(matrix.cols()),
                                                 static_cast<int64_t>(mat.nnz()),
                                                 mat.row_ptr(),
                                                 nullptr,
                                                 mat.col_ind(),
                                                 mat.val_ptr(),
                                                 CUDA_R_32I,
                                                 CUDA_P_TYPE,
                                                 CUDSS_MTYPE_GENERAL,
                                                 CUDSS_MVIEW_FULL,
                                                 CUDSS_BASE_ZERO),
                            "cudssMatrixCreateCsr");

        cuda::CudaVector rhs(n);
        cuda::CudaVector sol(n);
        cudssMatrix_t cudss_rhs {};
        cudssMatrix_t cudss_sol {};
        create_dense(rhs, sol, cudss_rhs, cudss_sol);

        runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, config, data, cudss_mat, cudss_sol, cudss_rhs),
                            "cudssExecute(analysis)");
        runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, config, data, cudss_mat, cudss_sol, cudss_rhs),
                            "cudssExecute(factorization)");

        destroy_dense(cudss_rhs, cudss_sol);
        timer.stop();
        factorization_time_ms = timer.elapsed();
    }

    ~CudssFactorizedSolver() {
        if (cudss_mat) runtime_check_cudss(cudssMatrixDestroy(cudss_mat), "cudssMatrixDestroy(mat)");
        if (data) runtime_check_cudss(cudssDataDestroy(handle, data), "cudssDataDestroy");
        if (config) runtime_check_cudss(cudssConfigDestroy(config), "cudssConfigDestroy");
        if (handle) runtime_check_cudss(cudssDestroy(handle), "cudssDestroy");
    }

    void solve(cuda::CudaVector& rhs, cuda::CudaVector& sol) {
        cudssMatrix_t cudss_rhs {};
        cudssMatrix_t cudss_sol {};
        create_dense(rhs, sol, cudss_rhs, cudss_sol);
        runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_SOLVE, config, data, cudss_mat, cudss_sol, cudss_rhs),
                            "cudssExecute(solve)");
        destroy_dense(cudss_rhs, cudss_sol);
    }

private:
    void create_dense(cuda::CudaVector& rhs,
                      cuda::CudaVector& sol,
                      cudssMatrix_t& cudss_rhs,
                      cudssMatrix_t& cudss_sol) {
        runtime_check_cudss(cudssMatrixCreateDn(&cudss_rhs,
                                                n,
                                                1,
                                                n,
                                                static_cast<CudaPrecision*>(rhs),
                                                CUDA_P_TYPE,
                                                CUDSS_LAYOUT_COL_MAJOR),
                            "cudssMatrixCreateDn(rhs)");
        runtime_check_cudss(cudssMatrixCreateDn(&cudss_sol,
                                                n,
                                                1,
                                                n,
                                                static_cast<CudaPrecision*>(sol),
                                                CUDA_P_TYPE,
                                                CUDSS_LAYOUT_COL_MAJOR),
                            "cudssMatrixCreateDn(sol)");
    }

    void destroy_dense(cudssMatrix_t& cudss_rhs, cudssMatrix_t& cudss_sol) {
        if (cudss_sol) {
            runtime_check_cudss(cudssMatrixDestroy(cudss_sol), "cudssMatrixDestroy(sol)");
            cudss_sol = {};
        }
        if (cudss_rhs) {
            runtime_check_cudss(cudssMatrixDestroy(cudss_rhs), "cudssMatrixDestroy(rhs)");
            cudss_rhs = {};
        }
    }
};

struct ForceBasisGpu {
    int n = 0;
    int n_basis = 0;
    cuda::CudaArray<CudaPrecision> alpha_values;
    std::vector<std::unique_ptr<cuda::CudaVector>> force_vectors;

    ForceBasisGpu(const NewmarkForceBasis& basis,
                  int n_in,
                  int n_steps,
                  double t_start,
                  double dt)
        : n(n_in),
          n_basis(static_cast<int>(basis.size())),
          alpha_values(static_cast<std::size_t>(basis.size()) * static_cast<std::size_t>(n_steps + 1))
    {
        std::vector<CudaPrecision> alpha_host(alpha_values.size(), CudaPrecision(0));
        force_vectors.reserve(basis.size());

        for (std::size_t basis_idx = 0; basis_idx < basis.size(); ++basis_idx) {
            const auto& amplitude = basis[basis_idx].first;
            const auto& force = basis[basis_idx].second;

            logging::error(force.size() == n, "Newmark GPU: force basis vector size mismatch.");

            auto force_gpu = std::make_unique<cuda::CudaVector>(n);
            force_gpu->upload(force.data());
            force_vectors.emplace_back(std::move(force_gpu));

            for (int step = 0; step <= n_steps; ++step) {
                const double time = t_start + step * dt;
                const Precision scale = amplitude ? amplitude->evaluate(time) : Precision(1);
                alpha_host[static_cast<std::size_t>(step) * basis.size() + basis_idx] =
                    static_cast<CudaPrecision>(scale);
            }
        }

        alpha_values.upload(alpha_host.data());
    }

    void assemble_step(int step, cuda::CudaVector& force) {
        force.clear();
        runtime_check_cuda(cublasSetPointerMode(cuda::manager.handle_cublas, CUBLAS_POINTER_MODE_DEVICE));
        for (int basis_idx = 0; basis_idx < n_basis; ++basis_idx) {
            const CudaPrecision* alpha =
                static_cast<CudaPrecision*>(alpha_values) + static_cast<std::size_t>(step) * n_basis + basis_idx;
            runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, alpha,
                                           *force_vectors[static_cast<std::size_t>(basis_idx)], 1,
                                           force, 1));
        }
        runtime_check_cuda(cublasSetPointerMode(cuda::manager.handle_cublas, CUBLAS_POINTER_MODE_HOST));
    }
};

void linear_combination(cuda::CudaVector& out,
                        cuda::CudaVector& x0,
                        CudaPrecision a0,
                        cuda::CudaVector& x1,
                        CudaPrecision a1,
                        cuda::CudaVector& x2,
                        CudaPrecision a2,
                        int n)
{
    out.copy(x0);
    runtime_check_cuda(CUBLAS_SCAL(cuda::manager.handle_cublas, n, &a0, out, 1));
    runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &a1, x1, 1, out, 1));
    runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &a2, x2, 1, out, 1));
}

void build_initial_acceleration(GpuSpmv& C_op,
                                GpuSpmv& K_op,
                                CudssFactorizedSolver* M_solver,
                                ForceBasisGpu& force_basis,
                                cuda::CudaVector& u,
                                cuda::CudaVector& v,
                                cuda::CudaVector& a,
                                cuda::CudaVector& rhs,
                                cuda::CudaVector& tmp,
                                const NewmarkIC& ic,
                                int n)
{
    if (ic.a0.size() == n) {
        a.upload(ic.a0.data());
        return;
    }

    logging::error(M_solver != nullptr, "Newmark GPU: M solver is required when initial acceleration is not supplied.");

    force_basis.assemble_step(0, rhs);
    C_op.apply(v, tmp, kOne, kZero);
    const CudaPrecision neg_one = CudaPrecision(-1);
    runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &neg_one, tmp, 1, rhs, 1));
    K_op.apply(u, tmp, kOne, kZero);
    runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &neg_one, tmp, 1, rhs, 1));
    M_solver->solve(rhs, a);
}

} // namespace

NewmarkResult
newmark_linear_gpu(const SparseMatrix& M,
                   const SparseMatrix& C,
                   const SparseMatrix& K,
                   const NewmarkIC& ic,
                   const NewmarkOpts& opts,
                   const NewmarkForceBasis& force_basis)
{
    logging::error(opts.dt > 0.0, "Newmark GPU: dt must be positive.");
    logging::error(opts.t_end >= opts.t_start, "Newmark GPU: t_end must be >= t_start.");
    logging::error(ic.u0.size() != 0 && ic.v0.size() != 0, "Newmark GPU: u0 and v0 must be provided.");
    logging::error(ic.u0.size() == ic.v0.size(), "Newmark GPU: u0 and v0 size mismatch.");
    logging::error(!force_basis.empty(), "Newmark GPU path requires force_basis.");

    const int n = checked_size(M.rows());
    logging::error(ic.u0.size() == n, "Newmark GPU: initial condition size mismatch.");

    cuda::manager.create_cuda();

    const double beta    = opts.beta;
    const double gamma   = opts.gamma;
    const double dt      = opts.dt;
    const double t_start = opts.t_start;
    const double t_end   = opts.t_end;

    const double a0 = 1.0 / (beta * dt * dt);
    const double a1 = gamma / (beta * dt);
    const double a2 = 1.0 / (beta * dt);
    const double a3 = 1.0 / (2.0 * beta) - 1.0;
    const double a4 = gamma / beta - 1.0;
    const double a5 = dt * (gamma / (2.0 * beta) - 1.0);

    const double Tspan = std::max(0.0, t_end - t_start);
    const int n_steps = static_cast<int>(std::ceil(Tspan / dt));
    const std::size_t n_history = static_cast<std::size_t>(n) * static_cast<std::size_t>(n_steps + 1);

    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "IMPLICIT NEWMARK-β GPU (fixed Δt) — linear transient");
    logging::info(true, "-----------------------------------------------------------------------------------------------");
    logging::info(true, "n dof  : ", n);
    logging::info(true, "basis  : ", static_cast<int>(force_basis.size()));
    logging::info(true, "nnz(M) : ", static_cast<long long>(M.nonZeros()));
    logging::info(true, "nnz(C) : ", static_cast<long long>(C.nonZeros()));
    logging::info(true, "nnz(K) : ", static_cast<long long>(K.nonZeros()));
    logging::info(true, "steps  : ", n_steps);
    logging::info(true, "dt      : ", std::scientific, std::setprecision(6), dt);
    logging::info(true, "t_start : ", std::scientific, std::setprecision(6), t_start);
    logging::info(true, "t_end   : ", std::scientific, std::setprecision(6), t_end);
    logging::info(true, "β, γ   : ", std::fixed, std::setprecision(4), beta, ", ", gamma);
    logging::up();

    SparseMatrix M_cpu = compressed_copy(M);
    SparseMatrix C_cpu = compressed_copy(C);
    SparseMatrix K_cpu = compressed_copy(K);

    logging::info(true, "");
    logging::info(true, "Building effective matrix A = K + a0*M + a1*C ...");
    Timer tBuild{};
    tBuild.start();
    SparseMatrix A_cpu = K_cpu + static_cast<Precision>(a0) * M_cpu + static_cast<Precision>(a1) * C_cpu;
    A_cpu.makeCompressed();
    tBuild.stop();
    logging::info(true, "A built: nnz(A) = ", static_cast<long long>(A_cpu.nonZeros()),
                       ", time = ", tBuild.elapsed(), " ms.");

    cuda::CudaVector u(n);
    cuda::CudaVector v(n);
    cuda::CudaVector a(n);
    cuda::CudaVector un(n);
    cuda::CudaVector vn(n);
    cuda::CudaVector an(n);
    cuda::CudaVector rhs(n);
    cuda::CudaVector tmp_m(n);
    cuda::CudaVector tmp_c(n);
    cuda::CudaVector tmp(n);

    u.upload(ic.u0.data());
    v.upload(ic.v0.data());

    ForceBasisGpu force_basis_gpu{force_basis, n, n_steps, t_start, dt};
    GpuSpmv M_op{M_cpu, tmp_m, rhs};
    GpuSpmv C_op{C_cpu, tmp_c, rhs};
    GpuSpmv K_op{K_cpu, u, tmp};

    logging::info(true, "Factorizing A with cuDSS (single factorization) ...");
    CudssFactorizedSolver A_solver{A_cpu, "A"};
    logging::info(true, "Factorization done in ", A_solver.factorization_time_ms, " ms.");

    std::unique_ptr<CudssFactorizedSolver> M_solver;
    if (ic.a0.size() != n) {
        logging::info(true, "Computing initial acceleration from equilibrium at t = ",
                             std::scientific, std::setprecision(6), t_start, " ...");
        logging::info(true, "Initial acceleration: factorizing mass matrix with cuDSS (one-time).");
        M_solver = std::make_unique<CudssFactorizedSolver>(M_cpu, "M");
    } else {
        logging::info(true, "Using user-supplied initial acceleration.");
    }

    Timer tInit{};
    tInit.start();
    build_initial_acceleration(C_op, K_op, M_solver.get(), force_basis_gpu, u, v, a, rhs, tmp, ic, n);
    tInit.stop();
    if (ic.a0.size() != n) {
        logging::info(true, "Initial acceleration: done in ", tInit.elapsed(), " ms.");
    }

    cuda::CudaArray<CudaPrecision> hist_u{n_history};
    cuda::CudaArray<CudaPrecision> hist_v{n_history};
    cuda::CudaArray<CudaPrecision> hist_a{n_history};

    copy_device_vector(hist_u, 0, u, static_cast<std::size_t>(n));
    copy_device_vector(hist_v, 0, v, static_cast<std::size_t>(n));
    copy_device_vector(hist_a, 0, a, static_cast<std::size_t>(n));

    const CudaPrecision ca0 = static_cast<CudaPrecision>(a0);
    const CudaPrecision ca1 = static_cast<CudaPrecision>(a1);
    const CudaPrecision ca2 = static_cast<CudaPrecision>(a2);
    const CudaPrecision ca3 = static_cast<CudaPrecision>(a3);
    const CudaPrecision ca4 = static_cast<CudaPrecision>(a4);
    const CudaPrecision ca5 = static_cast<CudaPrecision>(a5);
    const CudaPrecision cdt_one_minus_gamma = static_cast<CudaPrecision>(dt * (1.0 - gamma));
    const CudaPrecision cdt_gamma = static_cast<CudaPrecision>(dt * gamma);
    const CudaPrecision neg_one = CudaPrecision(-1);
    const CudaPrecision neg_a2 = -ca2;
    const CudaPrecision neg_a3 = -ca3;
    const int print_stride = std::max(1, n_steps / 20);

    logging::info(true, "");
    logging::info(true, "Starting Newmark GPU time marching (", n_steps, " steps) ...");

    Timer tAll{};
    tAll.start();
    Timer tLoop{};
    tLoop.start();

    for (int step = 0; step < n_steps; ++step) {
        const double tnp1 = t_start + (step + 1) * dt;

        force_basis_gpu.assemble_step(step + 1, rhs);

        linear_combination(tmp_m, u, ca0, v, ca2, a, ca3, n);
        linear_combination(tmp_c, u, ca1, v, ca4, a, ca5, n);

        M_op.apply(tmp_m, rhs, kOne, kOne);
        C_op.apply(tmp_c, rhs, kOne, kOne);

        A_solver.solve(rhs, un);

        an.copy(un);
        runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &neg_one, u, 1, an, 1));
        runtime_check_cuda(CUBLAS_SCAL(cuda::manager.handle_cublas, n, &ca0, an, 1));
        runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &neg_a2, v, 1, an, 1));
        runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &neg_a3, a, 1, an, 1));

        vn.copy(v);
        runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &cdt_one_minus_gamma, a, 1, vn, 1));
        runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, n, &cdt_gamma, an, 1, vn, 1));

        u.copy(un);
        v.copy(vn);
        a.copy(an);

        const std::size_t offset = static_cast<std::size_t>(step + 1) * static_cast<std::size_t>(n);
        copy_device_vector(hist_u, offset, u, static_cast<std::size_t>(n));
        copy_device_vector(hist_v, offset, v, static_cast<std::size_t>(n));
        copy_device_vector(hist_a, offset, a, static_cast<std::size_t>(n));

        if (((step + 1) % print_stride) == 0 || (step + 1) == n_steps) {
            logging::info(true, "Newmark step ",
                          std::setw(8), (step + 1), "/", n_steps,
                          "  t=", std::scientific, std::setprecision(6), tnp1, " s");
        }
    }

    tLoop.stop();
    tAll.stop();

    std::vector<CudaPrecision> hist_u_host(n_history);
    std::vector<CudaPrecision> hist_v_host(n_history);
    std::vector<CudaPrecision> hist_a_host(n_history);
    hist_u.download(hist_u_host.data());
    hist_v.download(hist_v_host.data());
    hist_a.download(hist_a_host.data());

    NewmarkResult out;
    out.t.reserve(n_steps + 1);
    out.u.reserve(n_steps + 1);
    out.v.reserve(n_steps + 1);
    out.a.reserve(n_steps + 1);

    for (int step = 0; step <= n_steps; ++step) {
        out.t.push_back(t_start + step * dt);

        DynamicVector u_step(n);
        DynamicVector v_step(n);
        DynamicVector a_step(n);
        const std::size_t offset = static_cast<std::size_t>(step) * static_cast<std::size_t>(n);
        for (int row = 0; row < n; ++row) {
            u_step[row] = static_cast<Precision>(hist_u_host[offset + row]);
            v_step[row] = static_cast<Precision>(hist_v_host[offset + row]);
            a_step[row] = static_cast<Precision>(hist_a_host[offset + row]);
        }

        out.u.push_back(std::move(u_step));
        out.v.push_back(std::move(v_step));
        out.a.push_back(std::move(a_step));
    }

    logging::info(true, "Time marching finished.");
    logging::up();
    logging::info(true, "Total wall time         : ", tAll.elapsed(), " ms");
    logging::info(true, "Loop time (steps only)  : ", tLoop.elapsed(), " ms");
    logging::info(true, "Avg per step            : ",
                         (n_steps > 0 ? tLoop.elapsed() / double(n_steps) : 0.0), " ms/step");
    logging::down();
    logging::down();

    return out;
}

#endif

} // namespace fem::solver::detail
