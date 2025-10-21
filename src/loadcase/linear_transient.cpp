/**
 * @file transient.cpp
 * @brief Linear transient analysis (implicit Newmark-β) using affine null-space map u = u_p + T q.
 */

#include "linear_transient.h"

#include "../core/logging.h"
#include "../core/timer.h"
#include "../reader/write_mtx.h"

#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"

#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"

#include <cmath>
#include <iomanip>
#include <algorithm>

namespace fem { namespace loadcase {

using fem::constraint::ConstraintTransformer;

Transient::Transient(ID id,
                     reader::Writer* writer_,
                     model::Model* model_)
    : LoadCase(id, writer_, model_) {}

void Transient::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR TRANSIENT ANALYSIS (Implicit Newmark-β)");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // (0) Materials/sections
    model->assign_sections();

    // (1) Unconstrained DOF index
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraint equations
    auto equations = Timer::measure(
        [&]() {
            auto groups = this->model->collect_constraints(active_dof_idx_mat, supps);
            //report_constraint_groups(groups); // enable if you want full print like buckling
            return groups.flatten();
        },
        "building constraints"
    );

    // (3) Global load matrix (node x 6) — time-independent for now
    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "building global load matrix"
    );

    // (4) Active stiffness K and mass M
    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );
    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M"
    );

    // (5) Reduce global loads -> active RHS f (n x 1)
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f"
    );

    // (6) Build constraint transformer (u = u_p + T q)
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions copt;
            copt.set.scale_columns = true; // robust QR, like buckling
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                K.rows(),            // n active DOFs
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
    logging::down();

    // (7) Reduce operators to master space: A = Tᵀ K T, Mr = Tᵀ M T
    auto A = Timer::measure(
        [&]() { return CT->assemble_A(K); },
        "assembling reduced stiffness A = T^T K T"
    );
    auto Mr = Timer::measure(
        [&]() { return CT->assemble_A(M); },
        "assembling reduced mass Mr = T^T M T"
    );

    // (8) Rayleigh damping in reduced space: Cr = α Mr + β A
    SparseMatrix Cr(A.rows(), A.cols());
    if (rayleigh.has_value()) {
        Cr = rayleigh->build(Mr, A);
    } else {
        Cr.setZero();
    }

    // (9) Reduced RHS (time-independent). With homogeneous supports u_p = 0:
    //     fr = Tᵀ (f - K u_p)  = Tᵀ f
    auto fr = Timer::measure(
        [&]() { return CT->assemble_b(K, f); }, // matches buckling style; with u_p=0 ⇒ Tᵀ f
        "reducing active RHS to master space (fr = T^T f)"
    );

    // Optional: dump matrices
    if (!stiffness_file.empty()) {
        write_mtx(stiffness_file + "_K.mtx", K , 0, 17);
        write_mtx(stiffness_file + "_A.mtx", A , 0, 17);
    }
    if (!mass_file.empty()) {
        write_mtx(mass_file + "_M.mtx", M  , 0, 17);
        write_mtx(mass_file + "_Mr.mtx", Mr, 0, 17);
    }
    if (!damping_file.empty()) {
        write_mtx(damping_file + "_C.mtx", Cr, 0, 17);
    }

    // (10) Build constant-force callback for the solver
    auto reduced_force = [fr](double /*t*/) -> DynamicVector { return fr; };

    // (11) Newmark options + IC (zeros in reduced space unless user-specified elsewhere)
    solver::NewmarkOpts optsNm;
    optsNm.beta   = beta;
    optsNm.gamma  = gamma;
    optsNm.dt     = dt;
    optsNm.t_end  = t_end;
    optsNm.device = device;
    optsNm.method = method;

    DynamicVector q0    = DynamicVector::Zero(A.rows());
    DynamicVector qdot0 = DynamicVector::Zero(A.rows());
    solver::NewmarkIC ic{q0, qdot0, DynamicVector()}; // a0 computed internally

    // (12) Solve reduced transient problem
    auto result = solver::newmark_linear(Mr, Cr, A, ic, optsNm, reduced_force);

    // (13) Write results with cadence
    const int n_steps = static_cast<int>(std::ceil(t_end / dt));
    int write_stride = std::max(1, write_every_steps);
    if (write_every_time > 0.0) {
        write_stride = std::max(1, static_cast<int>(std::round(write_every_time / dt)));
    }

    writer->add_loadcase(id);
    for (int k = 0; k <= n_steps; ++k) {
        if (k % write_stride != 0 && k != n_steps) continue; // always write last
        const auto& qk = result.u[static_cast<size_t>(k)];
        auto u_full = CT->recover_u(qk);
        auto U_mat  = mattools::expand_vec_to_mat(active_dof_idx_mat, u_full);
        writer->write_eigen_matrix(U_mat, "TRANSIENT_U_T" + std::to_string(k));
        // writer->write_eigen_matrix(result.v[k], "TRANSIENT_V_T"+std::to_string(k));
        // writer->write_eigen_matrix(result.a[k], "TRANSIENT_A_T"+std::to_string(k));
    }

    logging::info(true, "Transient analysis completed.");
}

}} // namespace fem::loadcase
