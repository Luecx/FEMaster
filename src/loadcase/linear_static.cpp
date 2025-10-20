/**
 * @file linear_static.cpp
 * @brief Implements the linear static load case leveraging constraint maps.
 *
 * The algorithm assembles the constrained system, reduces it via the null-space
 * map, solves for the reduced coordinates, and expands the solution and
 * reactions back to global DOFs.
 *
 * @see src/loadcase/linear_static.h
 * @see src/constraints/constraint_transformer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "linear_static.h"

#include "../constraints/builder/builder.h"
#include "../constraints/constraint_transformer.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../reader/write_mtx.h"
#include "../solve/eigen.h"

#include <utility>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

/**
 * @copydoc LinearStatic::LinearStatic
 */
LinearStatic::LinearStatic(ID id, reader::Writer* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

/**
 * @copydoc LinearStatic::run
 */
void LinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    model->assign_sections();

    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto equations = Timer::measure(
        [&]() {
            auto groups = model->collect_constraints(active_dof_idx_mat, supps);
            report_constraint_groups(groups);
            return groups.flatten();
        },
        "building constraints");

    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)");

    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K");

    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f");

    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions options;
            options.set.scale_columns = true;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                K.rows(),
                options);
        },
        "building constraint transformer");

    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)     : ", transformer->report().m);
    logging::info(true, "n (cols of C)     : ", transformer->report().n);
    logging::info(true, "rank(C)           : ", transformer->rank());
    logging::info(true, "masters (n-r)     : ", transformer->n_master());
    logging::info(true, "homogeneous       : ", transformer->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", transformer->feasible() ? "true" : "false");
    if (!transformer->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", transformer->report().residual_norm);
    }
    logging::down();

    auto A = Timer::measure(
        [&]() { return transformer->assemble_A(K); },
        "assembling reduced stiffness A = T^T K T");

    auto b = Timer::measure(
        [&]() { return transformer->assemble_b(K, f); },
        "assembling reduced RHS b = T^T (f - K u_p)");

    {
        bool bad_matrix = false;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<Precision>::InnerIterator it(A, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    bad_matrix = true;
                    break;
                }
            }
            if (bad_matrix) {
                break;
            }
        }
        logging::error(!bad_matrix, "Matrix A contains NaN/Inf entries");
        logging::error(b.allFinite(), "b contains NaN/Inf entries");
    }

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b); },
        "solving reduced system A q = b");

    auto u = Timer::measure(
        [&]() { return transformer->recover_u(q); },
        "recovering full displacement vector u");

    auto r = Timer::measure(
        [&]() { return transformer->reactions(K, f, q); },
        "computing reactions r = K u - f");

    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form");

    auto global_force_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r); },
        "expanding reactions to matrix form");

    NodeData stress;
    NodeData strain;
    std::tie(stress, strain) = Timer::measure(
        [&]() { return model->compute_stress_strain(global_disp_mat); },
        "interpolating stress and strain at nodes");

    if (!stiffness_file.empty()) {
        write_mtx(stiffness_file + "_K.mtx", K);
        write_mtx(stiffness_file + "_A.mtx", A);
        write_mtx_dense(stiffness_file + "_b.mtx", b);
    }

    writer->add_loadcase(id);
    writer->write_eigen_matrix(global_disp_mat, "DISPLACEMENT");
    writer->write_eigen_matrix(strain,          "STRAIN");
    writer->write_eigen_matrix(stress,          "STRESS");
    writer->write_eigen_matrix(global_load_mat, "DOF_LOADS");
    writer->write_eigen_matrix(global_force_mat,"NODAL_FORCES");

    transformer->post_check_static(K, f, u);
}

} // namespace loadcase
} // namespace fem
