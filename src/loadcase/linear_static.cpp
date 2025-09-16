/******************************************************************************
 * @file LinearStatic.cpp
 * @brief Linear static analysis using affine null-space constraints (u = u_p + T q).
 *
 * This implementation:
 *  - Builds the unconstrained system DOF indexing (active_dof_idx_mat).
 *  - Assembles global equations (constraints), global loads, and the active stiffness.
 *  - Builds the constraint transformer (C u = d  →  u = u_p + T q).
 *  - Reduces K and f: A = Tᵀ K T, b = Tᵀ (f - K u_p).
 *  - Solves A q = b, recovers u, computes reactions r = K u - f.
 *  - Expands u and r back to (nodes × 6) matrices for writing.
 *
 * Notes:
 *  - Loads are first reduced to the active DOFs via reduce_mat_to_vec(...) before use.
 *  - This path preserves SPD on A (for SPD K and feasible constraints).
 *  - The same transformer is re-usable across static/modal/buckling setups.
 *
 * @date Created on 04.09.2023
 ******************************************************************************/

#include "linear_static.h"

#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"

#include "../solve/eigen.h"
#include "../core/logging.h"

// New constraint stack
#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"
#include "../constraints/constraint_set.h"
#include "../constraints/constraint_builder.h"
#include "../constraints/constraint_map.h"

using fem::constraint::ConstraintTransformer;

fem::loadcase::LinearStatic::LinearStatic(ID id, reader::Writer* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

void fem::loadcase::LinearStatic::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC ANALYSIS");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    // (0) Section assignment (as before)
    m_model->assign_sections();

    // (1) Unconstrained system DOF index matrix (node×6 → active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Build constraint equations (supports, ties, couplings, …)
    auto equations = Timer::measure(
        [&]() { return this->m_model->build_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // (3) Global load matrix (node×6 layout); keep for output
    auto global_load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)"
    );

    // (4) Active global stiffness K (n×n), assembled using active_dof_idx_mat
    auto K = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );

    // (5) Reduce global loads to active RHS vector f (n×1)
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f"
    );

    // (6) Build constraint transformer (Set -> Builder -> Map)
    //     Uses the system DOF ids from the model and the active dimension n = K.rows().
    // Wrap creation of the transformer
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions copt;
            // Recommended: scale columns of C for robust QR:
            copt.set.scale_columns = true;
            // Optional tolerances:
            // copt.builder.rank_tol_rel = 1e-12;
            // copt.builder.feas_tol_rel = 1e-10;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,   // maps (node, dof) to global index
                K.rows(),             // total active DOFs n
                copt
            );
        },
        "building constraint transformer"
    );

    // Diagnostics
    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)     : ", CT->report().m);
    logging::info(true, "n (cols of C)     : ", CT->report().n);
    logging::info(true, "rank(C)           : ", CT->rank());
    logging::info(true, "masters (n-r)     : ", CT->n_master());
    logging::info(true, "homogeneous       : ", CT->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", CT->feasible() ? "true" : "false");
    if (!CT->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", CT->report().residual_norm);
    }
    logging::down();

    // (7) Assemble reduced system A q = b with A = Tᵀ K T, b = Tᵀ (f - K u_p)
    auto A = Timer::measure(
        [&]() { return CT->assemble_A(K); },
        "assembling reduced stiffness A = T^T K T"
    );
    auto b = Timer::measure(
        [&]() { return CT->assemble_b(K, f); },
        "assembling reduced RHS b = T^T (f - K u_p)"
    );

    // (8) Solve reduced system
    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b); },
        "solving reduced system A q = b"
    );

    // (9) Recover full displacement vector u (n×1)
    auto u = Timer::measure(
        [&]() { return CT->recover_u(q); },
        "recovering full displacement vector u"
    );

    // (10) Full reactions r = K u - f (n×1)
    auto r = Timer::measure(
        [&]() { return CT->reactions(K, f, q); },
        "computing reactions r = K u - f"
    );

    // (11) Expand u and r back to (node×6) for writer output
    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form"
    );
    auto global_force_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r); },
        "expanding reactions to matrix form"
    );

    // (12) Compute stresses and strains at the nodes (unchanged)
    NodeData stress, strain;
    std::tie(stress, strain) = Timer::measure(
        [&]() { return m_model->compute_stress_strain(global_disp_mat); },
        "Interpolating stress and strain at nodes"
    );

    // Optional: export K sub-block (active) if a path is configured
    if (!stiffness_file.empty()) {
        std::ofstream file(stiffness_file);
        for (int k = 0; k < K.outerSize(); ++k) {
            for (SparseMatrix::InnerIterator it(K, k); it; ++it) {
                file << it.row() << " " << it.col() << " " << it.value() << '\n';
            }
        }
        file.close();
    }

    // (13) Write results
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(global_disp_mat, "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain,          "STRAIN");
    m_writer->write_eigen_matrix(stress,          "STRESS");
    m_writer->write_eigen_matrix(global_load_mat, "DOF_LOADS");
    m_writer->write_eigen_matrix(global_force_mat,"NODAL_FORCES");

    // (14) Small consistency diagnostics (optional): projected residual
    {
        DynamicVector resid = K * u - f;
        DynamicVector red   = CT->map().apply_Tt(resid);
        logging::info(true, "");
        logging::info(true, "Post-checks");
        logging::up();
        logging::info(true, "||u||2                : ", u.norm());
        logging::info(true, "||C u - d||2          : ", (CT->set().C * u - CT->set().d).norm());
        logging::info(true, "||K u - f||2          : ", resid.norm());
        logging::info(true, "||T^T (K u - f)||2    : ", red.norm());
        logging::down();
    }
    logging::info(true, "");
    logging::info(true, "LINEAR STATIC ANALYSIS FINISHED");
    logging::info(true, "");
}
