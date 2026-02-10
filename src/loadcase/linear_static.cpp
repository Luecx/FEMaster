/**
 * @file linear_static.cpp
 * @brief Implements the linear static load case leveraging constraint maps.
 *
 * The algorithm assembles the constrained system, reduces it via the null-space
 * map, solves for the reduced coordinates, and expands the solution and
 * reactions back to global DOFs. Support reactions are recovered from
 * constraint multipliers (g_supp = C_supp^T λ_supp), ensuring nonzero values
 * even for supported DOFs without element stiffness (e.g. reference points
 * constrained via couplings).
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
#include "../mattools/reduce_mat_to_vec.h"
#include "../reader/write_mtx.h"
#include "../solve/eigen.h"
#include "../model/model.h"
#include "../model/element/element_structural.h"

#include "tools/inertia_relief.h"
#include "tools/rebalance_loads.h"

#include <utility>
#include <limits>
#include <algorithm>
#include <iomanip>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

/**
 * @copydoc LinearStatic::LinearStatic
 */
} // anonymous namespace

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

    // Build initial global load matrix (without inertia relief)
    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)");

    // ---------------------------------------------------------------------
    // Inertia relief: requires no supports in this load case
    // ---------------------------------------------------------------------
    if (inertia_relief) {
        rebalance_loads = true;
        logging::error(supps.empty(),
                       "InertiaRelief: cannot be used with *SUPPORT in this load case. "
                       "Remove all referenced support collectors.");

        Timer::measure(
            [&]() {
                fem::apply_inertia_relief(*model->_data, global_load_mat);

                // Add temporary RBM constraint (all nodes). Removed later after equations are built.
                model->add_rbm(std::string("NALL"));
            },
            "InertiaRelief: adjusting external load matrix and adding RBM");
    }

    if (rebalance_loads) {
        Timer::measure([&]() {fem::rebalance_loads(*model->_data, global_load_mat);},
            "rebalancing of loads");
    }

    // Build constraint groups (include RBM if inertia relief is enabled)
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

    // Remove temporary RBM again (equations already built)
    if (inertia_relief) {
        logging::error(model->_data != nullptr, "InertiaRelief: model data not initialized");
        logging::error(!model->_data->rbms.empty(),
                       "InertiaRelief: expected a temporary RBM to be present, but rbms is empty");
        model->_data->rbms.pop_back();
    }

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

    auto r_internal = Timer::measure(
        [&]() { return transformer->reactions(K, f, q); },
        "computing internal nodal forces r_int = K u - f");

    // Support reactions via constraint multipliers: g_supp = C_supp^T lambda_supp
    auto r_support = Timer::measure(
        [&]() { return transformer->support_reactions(K, f, q); },
        "computing support reactions via multipliers (C_supp^T lambda)");

    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form");

    auto global_react_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r_support); },
        "expanding support reactions to matrix form");

    auto section_force_mat = Timer::measure(
        [&]() { return model->compute_section_forces(global_disp_mat); },
        "computing beam section forces");

    auto [stress, strain] = Timer::measure(
        [&]() { return model->compute_stress_strain(global_disp_mat); },
        "interpolating stress and strain at nodes");

    if (!stiffness_file.empty()) {
        write_mtx(stiffness_file + "_K.mtx", K);
        write_mtx(stiffness_file + "_A.mtx", A);
        write_mtx_dense(stiffness_file + "_b.mtx", b);
    }

    // Build a support-only mask for reaction reporting: NaN everywhere except
    // DOFs that were constrained by supports.
    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);
    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 && e.node_id < support_mask.rows() &&
                e.dof >= 0 && e.dof < support_mask.cols())
            {
                // Only mark if DOF exists in the active numbering
                if (active_dof_idx_mat(e.node_id, e.dof) != -1) {
                    support_mask(e.node_id, e.dof) = true;
                }
            }
        }
    }

    model::Field reaction_masked{"REACTION_FORCES", model::FieldDomain::NODE,
                                 global_react_mat.rows,
                                 global_react_mat.components};
    reaction_masked.fill_nan();
    for (int i = 0; i < reaction_masked.rows; ++i) {
        for (int j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = global_react_mat(i, j);
            }
        }
    }

    writer->add_loadcase(id);
    writer->write_field(global_disp_mat         , "DISPLACEMENT");
    writer->write_field(strain                  , "STRAIN");
    writer->write_field(stress                  , "STRESS");
    writer->write_field(global_load_mat         , "EXTERNAL_FORCES");
    writer->write_field(reaction_masked         , "REACTION_FORCES");
    writer->write_eigen_matrix(section_force_mat, "LOCAL_SECTION_FORCES", 2);

    transformer->post_check_static(K, f, u);
}

} // namespace loadcase
} // namespace fem