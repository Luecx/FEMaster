/**
 * @file transient.cpp
 * @brief Linear transient analysis (implicit Newmark-β) using affine null-space map u = u_p + T q.
 */

#include "linear_transient.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/equation.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../writer/write_mtx.h"

#include <algorithm>
#include <cmath>
#include <iomanip>

namespace fem { namespace loadcase {

using fem::constraint::ConstraintTransformer;

Transient::Transient(ID id,
                     io::writer::ResultWriters* writer_,
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
    model->step_begin();

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

    // (3) Active stiffness K and mass M
    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K"
    );
    auto M = Timer::measure(
        [&]() { return model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix M"
    );

    // (4) Build constraint transformer (u = u_p + T q)
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::Options copt;
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
    logging::info(true, "m (rows of C)     : ", CT->report().equations);
    logging::info(true, "n (cols of C)     : ", CT->report().dofs);
    logging::info(true, "rank(C)           : ", CT->rank());
    logging::info(true, "solver unknowns    : ", CT->unknowns());
    logging::info(true, "homogeneous       : ", CT->homogeneous() ? "true" : "false");
    logging::down();

    // (5) Reduce operators to master space: A = Tᵀ K T, Mr = Tᵀ M T
    auto A = Timer::measure(
        [&]() { return CT->assemble_system_matrix(K); },
        "assembling reduced stiffness A = T^T K T"
    );
    auto Mr = Timer::measure(
        [&]() { return CT->assemble_system_matrix(M); },
        "assembling reduced mass Mr = T^T M T"
    );

    // (6) Rayleigh damping in reduced space: Cr = α Mr + β A
    SparseMatrix Cr(A.rows(), A.cols());
    if (rayleigh.has_value()) {
        Cr = rayleigh->build(Mr, A);
    } else {
        Cr.setZero();
    }

    // Optional: dump matrices
    if (!stiffness_file.empty()) {
        io::writer::write_mtx(stiffness_file + "_K.mtx", K , 0, 17);
        io::writer::write_mtx(stiffness_file + "_A.mtx", A , 0, 17);
    }
    if (!mass_file.empty()) {
        io::writer::write_mtx(mass_file + "_M.mtx", M  , 0, 17);
        io::writer::write_mtx(mass_file + "_Mr.mtx", Mr, 0, 17);
    }
    if (!damping_file.empty()) {
        io::writer::write_mtx(damping_file + "_C.mtx", Cr, 0, 17);
    }

    // (7) Build time-dependent reduced force callback for the solver
    auto reduced_force = [this,
                           &active_dof_idx_mat,
                           CT_ptr = CT.get(),
                           &K](double time) -> DynamicVector {
        auto load_matrix = model->build_load_matrix(this->loads, time);
        auto f_active = mattools::reduce_mat_to_vec(active_dof_idx_mat, load_matrix);
        return CT_ptr->assemble_system_rhs(K, f_active);
    };

    solver::NewmarkForceBasis reduced_force_basis;
    if (device == solver::GPU) {
        auto load_basis = model->build_load_basis(this->loads);
        reduced_force_basis.reserve(load_basis.size());
        for (auto& [amplitude, load_matrix] : load_basis) {
            auto f_active = mattools::reduce_mat_to_vec(active_dof_idx_mat, load_matrix);
            reduced_force_basis.emplace_back(amplitude, CT->assemble_system_rhs(K, f_active));
        }
    }

    // (8) Newmark options + IC (zeros in reduced space unless user provided an initial velocity field)
    solver::NewmarkOpts optsNm;
    optsNm.beta    = beta;
    optsNm.gamma   = gamma;
    optsNm.dt      = dt;
    optsNm.t_start = t_start;
    optsNm.t_end   = t_end;
    optsNm.device  = device;
    optsNm.method  = method;

    DynamicVector q0    = DynamicVector::Zero(A.rows());
    DynamicVector qdot0 = DynamicVector::Zero(A.rows());

    // Optional: initial velocity from node FIELD (must have exactly 6 components)
    if (initial_velocity) {
        logging::error(initial_velocity->domain == fem::model::FieldDomain::NODE, "Initial velocity field must be a node field");
        logging::error(initial_velocity->components == 6, "Initial velocity field must have exactly 6 components");

        auto v_active = mattools::reduce_mat_to_vec(active_dof_idx_mat, *initial_velocity);
        // Map to reduced initial velocity qdot0 = T^T v_active
        qdot0 = CT->project_vector(v_active);
    }
    solver::NewmarkIC ic{q0, qdot0, DynamicVector()}; // a0 computed internally

    // (9) Solve reduced transient problem
    auto result = solver::newmark_linear(Mr, Cr, A, ic, optsNm, reduced_force, reduced_force_basis);

    // (10) Write results with cadence
    const int n_steps = static_cast<int>(std::ceil(std::max(0.0, t_end - t_start) / dt));
    int write_stride = std::max(1, write_every_steps);
    if (write_every_time > 0.0) {
        write_stride = std::max(1, static_cast<int>(std::round(write_every_time / dt)));
    }

    writer->add_loadcase(id, io::writer::WriterStepType::Dynamic);
    for (int k = 0; k <= n_steps; ++k) {
        if (k % write_stride != 0 && k != n_steps) continue; // always write last
        const auto& qk  = result.u[static_cast<size_t>(k)];
        const auto& qvk = result.v[static_cast<size_t>(k)];
        const auto& qak = result.a[static_cast<size_t>(k)];
        const Precision frame_time =
            static_cast<Precision>(t_start + static_cast<double>(k) * dt);

        // Displacement: u = T q + u_p
        auto u_full = CT->recover_displacement(qk);
        auto U_mat  = mattools::expand_vec_to_mat(active_dof_idx_mat, u_full);
        writer->write_field(U_mat, "DISPLACEMENT_" + std::to_string(k), model->_data.get(), frame_time);

        // Velocity recovery uses the active constraint backend mapping.
        auto v_full = CT->recover_velocity(qvk);
        auto V_mat  = mattools::expand_vec_to_mat(active_dof_idx_mat, v_full);
        writer->write_field(V_mat, "VELOCITY_" + std::to_string(k), model->_data.get(), frame_time);

        // Acceleration recovery uses the active constraint backend mapping.
        auto a_full = CT->recover_acceleration(qak);
        auto A_mat  = mattools::expand_vec_to_mat(active_dof_idx_mat, a_full);
        writer->write_field(A_mat, "ACCELERATION_" + std::to_string(k), model->_data.get(), frame_time);
    }

    logging::info(true, "Transient analysis completed.");
    model->step_end();
}
}} // namespace fem::loadcase
