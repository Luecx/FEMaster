#include "solve_sparse_indirect.h"

#include "../../core/logging.h"

#ifdef SUPPORT_GPU
#ifdef USE_AMGX
#include <amgx_c.h>

#include <Eigen/SparseCore>

#include <chrono>
#include <limits>
#include <string>
#endif
#endif

namespace fem::solver::detail {

#ifdef SUPPORT_GPU
#ifdef USE_AMGX
namespace {

constexpr const char* AMGX_CONFIG =
    "config_version=2,"
    "determinism_flag=1,"
    "solver(main)=FGMRES,"
    "main:scaling=DIAGONAL_SYMMETRIC,"
    "main:preconditioner(amg)=AMG,"
    "main:gmres_n_restart=32,"
    "main:use_scalar_norm=1,"
    "amg:algorithm=AGGREGATION,"
    "amg:selector=SIZE_2,"
    "amg:smoother(smooth)=MULTICOLOR_DILU,"
    "smooth:relaxation_factor=0.75,"
    "amg:max_uncolored_percentage=0.05,"
    "amg:matrix_coloring_scheme=PARALLEL_GREEDY,"
    "amg:presweeps=0,"
    "amg:postsweeps=3,"
    "amg:interpolator=D2,"
    "amg:coarse_solver=DENSE_LU_SOLVER,"
    "amg:min_coarse_rows=32,"
    "amg:max_iters=1,"
    "amg:max_levels=50,"
    "amg:cycle=V,"
    "main:max_iters=1000,"
    "main:monitor_residual=1,"
    "main:store_res_history=1,"
    "main:convergence=RELATIVE_INI,"
    "main:tolerance=1e-10,"
    "main:norm=L2,"
    "main:print_solve_stats=0,"
    "main:obtain_timings=1,"
    "amg:print_grid_stats=0";

void check_amgx(AMGX_RC rc, const char* call) {
    if (rc == AMGX_RC_OK) {
        return;
    }

    char message[4096] = {};
    AMGX_get_error_string(rc, message, sizeof(message));
    logging::error(false, call, " failed: ", message);
}

class AmgxRuntime {
public:
    AmgxRuntime() {
        check_amgx(AMGX_initialize(), "AMGX_initialize");
    }

    ~AmgxRuntime() {
        AMGX_finalize();
    }

    AmgxRuntime(const AmgxRuntime&) = delete;
    AmgxRuntime& operator=(const AmgxRuntime&) = delete;
};

struct AmgxHandles {
    AMGX_config_handle config = nullptr;
    AMGX_resources_handle resources = nullptr;
    AMGX_matrix_handle matrix = nullptr;
    AMGX_vector_handle rhs = nullptr;
    AMGX_vector_handle solution = nullptr;
    AMGX_solver_handle solver = nullptr;

    ~AmgxHandles() {
        if (solver) {
            AMGX_solver_destroy(solver);
        }
        if (solution) {
            AMGX_vector_destroy(solution);
        }
        if (rhs) {
            AMGX_vector_destroy(rhs);
        }
        if (matrix) {
            AMGX_matrix_destroy(matrix);
        }
        if (resources) {
            AMGX_resources_destroy(resources);
        }
        if (config) {
            AMGX_config_destroy(config);
        }
    }
};

AmgxRuntime& amgx_runtime() {
    static AmgxRuntime runtime;
    return runtime;
}

const char* status_name(AMGX_SOLVE_STATUS status) {
    switch (status) {
        case AMGX_SOLVE_SUCCESS:       return "success";
        case AMGX_SOLVE_FAILED:        return "failed";
        case AMGX_SOLVE_DIVERGED:      return "diverged";
        case AMGX_SOLVE_NOT_CONVERGED: return "not converged";
    }
    return "unknown";
}

double elapsed_ms(const std::chrono::steady_clock::time_point& begin) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
}

void log_residual_history(AMGX_solver_handle solver, int iterations) {
    for (int iteration = 0; iteration < iterations; ++iteration) {
        if (iteration >= 10 && (iteration + 1) % 100 != 0 && iteration != iterations - 1) {
            continue;
        }

        double residual = 0.0;
        check_amgx(AMGX_solver_get_iteration_residual(solver, iteration, 0, &residual),
                   "AMGX_solver_get_iteration_residual");
        logging::info(true, "AMGX iteration ", iteration + 1,
                      " residual: ", residual);
    }
}

} // namespace

DynamicMatrix solve_indirect_gpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs_matrix) {
    (void)amgx_runtime();

    logging::error(mat.rows() <= std::numeric_limits<int>::max(),
                   "AMGX currently requires 32-bit row indices");
    logging::error(mat.nonZeros() <= std::numeric_limits<int>::max(),
                   "AMGX currently requires 32-bit nonzero indices");

    using RowSparseMatrix = Eigen::SparseMatrix<Precision, Eigen::RowMajor, int>;
    RowSparseMatrix csr = mat;
    csr.makeCompressed();

    const int n = static_cast<int>(csr.rows());
    const int nnz = static_cast<int>(csr.nonZeros());

#ifdef DOUBLE_PRECISION
    constexpr AMGX_Mode mode = AMGX_mode_dDDI;
#else
    constexpr AMGX_Mode mode = AMGX_mode_dFFI;
#endif

    logging::info(true, "AMGX backend: FGMRES + aggregation AMG + multicolor DILU + symmetric diagonal scaling");
    logging::info(true, "AMGX upload: N=", n, ", nnz=", nnz);

    DynamicMatrix result(mat.rows(), rhs_matrix.cols());
    double upload_ms = 0.0;
    double setup_ms = 0.0;
    double solve_ms = 0.0;
    Eigen::Index failed_column = -1;
    AMGX_SOLVE_STATUS failed_status = AMGX_SOLVE_SUCCESS;

    {
        AmgxHandles handles;
        check_amgx(AMGX_config_create(&handles.config, AMGX_CONFIG),
                   "AMGX_config_create");
        check_amgx(AMGX_resources_create_simple(&handles.resources, handles.config),
                   "AMGX_resources_create_simple");
        check_amgx(AMGX_matrix_create(&handles.matrix, handles.resources, mode),
                   "AMGX_matrix_create");
        check_amgx(AMGX_vector_create(&handles.rhs, handles.resources, mode),
                   "AMGX_vector_create(rhs)");
        check_amgx(AMGX_vector_create(&handles.solution, handles.resources, mode),
                   "AMGX_vector_create(solution)");
        check_amgx(AMGX_solver_create(&handles.solver, handles.resources, mode, handles.config),
                   "AMGX_solver_create");

        const auto upload_begin = std::chrono::steady_clock::now();
        check_amgx(AMGX_matrix_upload_all(handles.matrix,
                                          n,
                                          nnz,
                                          1,
                                          1,
                                          csr.outerIndexPtr(),
                                          csr.innerIndexPtr(),
                                          csr.valuePtr(),
                                          nullptr),
                   "AMGX_matrix_upload_all");
        check_amgx(AMGX_vector_bind(handles.rhs, handles.matrix),
                   "AMGX_vector_bind(rhs)");
        check_amgx(AMGX_vector_bind(handles.solution, handles.matrix),
                   "AMGX_vector_bind(solution)");
        upload_ms = elapsed_ms(upload_begin);

        const auto setup_begin = std::chrono::steady_clock::now();
        check_amgx(AMGX_solver_setup(handles.solver, handles.matrix),
                   "AMGX_solver_setup");
        setup_ms = elapsed_ms(setup_begin);

        for (Eigen::Index column = 0; column < rhs_matrix.cols(); ++column) {
            DynamicVector rhs = rhs_matrix.col(column);
            DynamicVector solution = DynamicVector::Zero(mat.rows());

            check_amgx(AMGX_vector_upload(handles.rhs, n, 1, rhs.data()),
                       "AMGX_vector_upload(rhs)");
            check_amgx(AMGX_vector_set_zero(handles.solution, n, 1),
                       "AMGX_vector_set_zero(solution)");

            const auto solve_begin = std::chrono::steady_clock::now();
            check_amgx(AMGX_solver_solve(handles.solver, handles.rhs, handles.solution),
                       "AMGX_solver_solve");
            solve_ms += elapsed_ms(solve_begin);

            AMGX_SOLVE_STATUS status = AMGX_SOLVE_FAILED;
            int iterations = 0;
            check_amgx(AMGX_solver_get_status(handles.solver, &status),
                       "AMGX_solver_get_status");
            check_amgx(AMGX_solver_get_iterations_number(handles.solver, &iterations),
                       "AMGX_solver_get_iterations_number");

            log_residual_history(handles.solver, iterations);
            logging::info(true, "RHS ", column, ": ", iterations,
                          " AMGX iterations, status=", status_name(status));

            if (status != AMGX_SOLVE_SUCCESS) {
                failed_column = column;
                failed_status = status;
                break;
            }

            check_amgx(AMGX_vector_download(handles.solution, solution.data()),
                       "AMGX_vector_download(solution)");
            result.col(column) = solution;
        }
    }

    logging::info(true, "AMGX matrix upload: ", upload_ms, " ms");
    logging::info(true, "AMGX AMG setup    : ", setup_ms, " ms");
    logging::info(true, "AMGX solve total  : ", solve_ms, " ms");

    logging::error(failed_column < 0,
                   "AMGX solve for RHS ", failed_column, " ", status_name(failed_status));

    return result;
}

#else

DynamicMatrix solve_indirect_gpu(SparseMatrix&,
                                 const DynamicMatrix&) {
    logging::error(false,
                   "GPU indirect solving requires AMGX integration on this branch");
    return {};
}

#endif
#else

DynamicMatrix solve_indirect_gpu(SparseMatrix&,
                                 const DynamicMatrix&) {
    logging::error(false, "GPU indirect solving is not available in this build");
    return {};
}

#endif

} // namespace fem::solver::detail
