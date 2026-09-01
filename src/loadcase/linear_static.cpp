/**
 * @file linear_static.cpp
 * @brief Implements the linear static load case leveraging constraint maps.
 */

#include "linear_static.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/rbm.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/element/element_structural.h"
#include "../model/model.h"
#include "../solve/eigval/solve_eigval.h"
#include "../io/writer/write_mtx.h"
#include "tools/inertia_relief.h"
#include "tools/rebalance_loads.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

void LinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    model->assign_sections();
    model->step_begin();

    // DEBUG A/B: force the current nodal positions back to the reference geometry
    // before the linear stiffness matrix is assembled. The old solid GL path then
    // evaluates exactly at F = I if POSITION was carrying any prior/current state.
    *model->_data->positions = *model->_data->positions_reference;

    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)");

    if (inertia_relief) {
        logging::error(supps.empty(),
            "InertiaRelief: cannot be used with *SUPPORT in this load case. "
            "Remove all referenced support collectors.");

        Timer::measure(
            [&]() {
                fem::apply_inertia_relief(
                    *model->_data,
                    global_load_mat,
                    inertia_relief_consider_point_masses
                );

                logging::error(model->_data->elem_sets.has("EALL")
                            && model->_data->elem_sets.get("EALL") != nullptr,
                    "InertiaRelief: EALL element set is not available");

                // Inertia relief adds one temporary RBM directly to the same
                // collection used by parser-created constraints. It is removed
                // again immediately after constraint collection below.
                model->_data->rbms.emplace_back(model->_data->elem_sets.get("EALL"));
            },
            "InertiaRelief: adjusting external load matrix and adding RBM");
    }

    if (rebalance_loads) {
        logging::error(supps.empty(),
            "Rebalancing Loads: cannot be used with *SUPPORT in this load case. "
            "Remove all referenced support collectors.");

        Timer::measure(
            [&]() { fem::rebalance_loads(*model->_data, global_load_mat); },
            "rebalancing of loads");
    }

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints");

    report_constraint_groups(groups);
    auto equations = groups.flatten();

    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K");

    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f");

    if (constraint_method == ConstraintTransformer::Method::Lagrange && method == solver::INDIRECT) {
        logging::error(false,
            "Invalid solver/constraint combination\n"
            "Constraint | Backend   | DIRECT       | INDIRECT\n"
            "NULLSPACE  | CPU MKL   | Yes          | Yes\n"
            "NULLSPACE  | CPU Eigen | Yes          | Yes\n"
            "NULLSPACE  | GPU       | Yes          | Yes\n"
            "NULLSPACE  | GPU cuDSS | Yes          | Yes\n"
            "LAGRANGE   | CPU MKL   | Yes          | No\n"
            "LAGRANGE   | CPU Eigen | Limited      | No\n"
            "LAGRANGE   | GPU       | No           | No\n"
            "LAGRANGE   | GPU cuDSS | Yes          | No\n"
            "ELIMINATION| CPU MKL   | Yes          | Yes\n"
            "ELIMINATION| CPU Eigen | Yes          | Yes\n"
            "ELIMINATION| GPU       | Yes          | Yes\n"
            "ELIMINATION| GPU cuDSS | Yes          | Yes");
    }

    const auto direct_matrix_type =
        constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;
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
    logging::info(true, "m (rows of C)     : ", transformer->report().equations);
    logging::info(true, "n (cols of C)     : ", transformer->report().dofs);
    if (transformer->rank_known()) {
        logging::info(true, "rank(C)           : ", transformer->rank());
    } else {
        logging::info(true, "rank(C)           : not computed");
    }
    logging::info(true, "method            : ", transformer->method_name());
    logging::info(true, "solver unknowns   : ", transformer->unknowns());
    logging::info(true, "homogeneous       : ", transformer->homogeneous() ? "true" : "false");
    logging::info(true, "feasible          : ", transformer->feasible() ? "true" : "false");
    if (!transformer->feasible()) {
        logging::info(true, "residual ||C u - d|| : ", transformer->report().residual_norm);
    }
    logging::down();

    if (inertia_relief) {
        logging::error(!model->_data->rbms.empty(),
            "InertiaRelief: expected a temporary RBM to be present, but rbms is empty");
        model->_data->rbms.pop_back();
    }

    auto A = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(K); },
        "assembling constraint system matrix");

    auto b = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(K, f); },
        "assembling constraint system RHS");

    {
        bool bad_matrix = false;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<Precision>::InnerIterator it(A, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    bad_matrix = true;
                    break;
                }
            }
            if (bad_matrix) break;
        }
        logging::error(!bad_matrix,
            "Matrix A contains NaN/Inf entries");
        logging::error(b.allFinite(),
            "b contains NaN/Inf entries");
    }

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b, direct_matrix_type); },
        "solving constraint system");

    auto u = Timer::measure(
        [&]() { return transformer->recover_displacement(q); },
        "recovering full displacement vector u");

    auto r_internal = Timer::measure(
        [&]() { return transformer->reactions(K, f, q); },
        "computing internal nodal forces r_int = K u - f");
    (void) r_internal;

    auto r_support = Timer::measure(
        [&]() { return transformer->support_reactions(K, f, q); },
        "computing support reactions via multipliers (C_supp^T lambda)");

    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form");

    auto global_react_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r_support); },
        "expanding support reactions to matrix form");

    auto section_forces = Timer::measure(
        [&]() { return model->compute_section_forces(global_disp_mat); },
        "computing beam section forces");

    auto shear_flow = Timer::measure(
        [&]() { return model->compute_shear_flow(global_disp_mat); },
        "computing shear-flow output");

    auto stress_strain = Timer::measure(
        [&]() { return model->compute_stress_nodal(global_disp_mat, false); },
        "interpolating stress and strain at nodes");
    auto stress = std::move(std::get<0>(stress_strain));
    auto strain = std::move(std::get<1>(stress_strain));

    auto stress_top_bottom = Timer::measure(
        [&]() { return model->compute_stress_top_bot(global_disp_mat, false); },
        "interpolating top/bottom stress at nodes");
    auto stress_top = std::move(std::get<0>(stress_top_bottom));
    auto stress_bot = std::move(std::get<1>(stress_top_bottom));

    auto shell_resultants = Timer::measure(
        [&]() { return model->compute_shell_resultants(global_disp_mat); },
        "interpolating shell resultants at nodes");

    if (!stiffness_file.empty()) {
        io::writer::write_mtx(stiffness_file + "_K.mtx", K);
        io::writer::write_mtx(stiffness_file + "_A.mtx", A);
        io::writer::write_mtx_dense(stiffness_file + "_b.mtx", b);
    }

    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);
    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 && e.node_id < support_mask.rows()
             && e.dof < support_mask.cols()
             && active_dof_idx_mat(e.node_id, e.dof) != -1) {
                support_mask(e.node_id, e.dof) = true;
            }
        }
    }

    model::Field reaction_masked{
        "REACTION_FORCES",
        model::FieldDomain::NODE,
        global_react_mat.rows,
        global_react_mat.components
    };
    reaction_masked.fill_nan();
    for (Index i = 0; i < reaction_masked.rows; ++i) {
        for (Index j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) reaction_masked(i, j) = global_react_mat(i, j);
        }
    }

    Timer::measure(
        [&]() {
            writer->add_loadcase(id, io::writer::WriterStepType::Static);
            writer->write_field(global_disp_mat , "DISPLACEMENT", model->_data.get());
            writer->write_field(strain          , "STRAIN", model->_data.get());
            writer->write_field(stress          , "STRESS", model->_data.get());
            writer->write_field(stress_top      , "STRESS_TOP", model->_data.get());
            writer->write_field(stress_bot      , "STRESS_BOT", model->_data.get());
            writer->write_field(shell_resultants, "SHELL_RESULTANTS", model->_data.get());
            writer->write_field(global_load_mat , "EXTERNAL_FORCES", model->_data.get());
            writer->write_field(reaction_masked , "REACTION_FORCES", model->_data.get());
            writer->write_field(section_forces  , "LOCAL_SECTION_FORCES", model->_data.get());
            if (shear_flow.rows > 0) {
                writer->write_field(shear_flow, "SHEAR_FLOW", model->_data.get());
            }
        },
        "writing result fields");

    transformer->post_check_static(K, f, q);
    model->step_end();
}

} // namespace loadcase
} // namespace fem
