/**
 * @file LinearStaticTopo.cpp
 * @brief Linear static analysis with topology optimization (density/orientation).
 *
 * Overview
 * --------
 * This load case performs a linear static solve while exposing per-element
 * topology parameters (density and orientation) that scale the stiffness
 * assembly. Constraints are enforced via an affine null-space parameterization
 *   u = u_p + T q
 * so that the reduced system
 *   A q = b,  with  A = T^T K T,  b = T^T (f - K u_p)
 * is unconstrained and (for SPD K and feasible constraints) symmetric positive
 * definite. This matches the LinearStatic pipeline but inserts the SIMP-like
 * stiffness scaling and orientation before assembling K.
 *
 * Pipeline (null-space constrained, similar to LinearStatic)
 * ----------------------------------------------------------
 *  (0) Assign sections/materials; push SIMP-like element scaling (rho^p) and
 *      per-element orientation angles to the model data store.
 *  (1) Build unconstrained DOF indexing (node x 6 -> active dof id or -1).
 *  (2) Build constraint equations (supports/ties/couplings -> C u = d).
 *  (3) Assemble global loads (node x 6), reduce loads to active RHS f.
 *  (4) Assemble active stiffness K(rho, p, theta) using current density/orientation.
 *  (5) Build ConstraintTransformer -> u = u_p + T q (wrapped for timing).
 *  (6) Reduce and solve: A = T^T K T, b = T^T (f - K u_p); solve A q = b.
 *  (7) Recover u and support reactions via multipliers (g_supp = C_supp^T λ_supp);
 *      expand to node x 6.
 *  (8) Post-processing: nodal stress/strain.
 *  (9) Topology metrics: compliance and sensitivities (density/orientation), volume.
 * (10) Write results.
 * (11) Cleanup temporary element data from the model store.
 * (12) Optional diagnostics: projected residual checks.
 *
 * Notes
 * -----
 * - Density scaling uses a SIMP-like exponent "exponent" (rho^p). For p > 1,
 *   intermediate densities are penalized. Orientation is passed as angles per
 *   element and consumed by the model during K assembly.
 * - Because constraints are handled via T, the reduced matrix A remains SPD
 *   when K is SPD and the constraint set is feasible (no contradictory rows).
 * - Direct solvers benefit from assembling A explicitly. If memory becomes a
 *   concern on very large problems, consider operator mode (matrix-free).
 *
 * @date    27.08.2024
 * @author  Finn Eggers
 */

#include "linear_static_topo.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/equation.h"
#include "../core/logging.h"
#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_vec.h"
#include "../solve/eigval/solve_eigval.h"
#include "tools/inertia_relief.h"
#include "tools/rebalance_loads.h"

#include <limits>
#include <cmath>

using fem::constraint::ConstraintTransformer;

namespace fem { namespace loadcase {

/**
 * @class LinearStaticTopo
 * @brief Linear static analysis with topology parameters (density/orientation).
 *
 * Members
 * -------
 * - density    : per-element densities (initialized to 1.0).
 * - orientation: per-element orientation angles (3 per element, initialized to 0).
 */
LinearStaticTopo::LinearStaticTopo(ID id, io::writer::ResultWriters* writer, model::Model* model)
    : LinearStatic(id, writer, model) {}

/**
 * @brief Execute the topology-aware linear static analysis.
 *
 * Implementation notes
 * --------------------
 * - We mirror the logging/timing style from LinearStatic: all heavyweight
 *   stages are wrapped with Timer::measure(...) and produce consistent labels.
 * - ConstraintTransformer creation is wrapped for timing; diagnostics are
 *   printed (rank, homogeneity, feasibility).
 */
void LinearStaticTopo::run() {
    // Banner
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC TOPO");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // (0) Sections/materials
    model->assign_sections();

    // Inject topology parameters into the model's element data store so the
    // stiffness builder can read them during assembly.
    logging::error(density      != nullptr, "LinearStaticTopo: density field not initialized");

    // copy the density to the stiffness field
    model::Field::Ptr stiffness = std::make_shared<model::Field>(
        "TOPO_DENSITY_STIFFNESS", model::FieldDomain::ELEMENT, density->rows, 1);

    for (Index i = 0; i < stiffness->rows; i++) {
        if (!std::isnan((*density)(i)) && std::isfinite((*density)(i)))
            (*stiffness)(i) = std::pow((*density)(i), this->exponent);
    }

    // assign the model stiffness and
    model->_data->element_stiffness_scale = stiffness;
    model->_data->material_orientation    = orientation;
    model->step_begin();

    // (1) Unconstrained DOF indexing (node x 6 -> active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Global loads (node x 6)
    auto global_load_mat = Timer::measure(
        [&]() { return this->model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)"
    );

    // (2a) Inertia relief (optional) and temporary RBM
    if (inertia_relief) {
        logging::error(supps.empty(),
                       "InertiaRelief: cannot be used with *SUPPORT in this load case. "
                       "Remove all referenced support collectors.");

        Timer::measure(
            [&]() {
                fem::apply_inertia_relief(*model->_data,
                                          global_load_mat,
                                          inertia_relief_consider_point_masses);
                // Add temporary RBM constraint (all elements). Removed later after equations are built.
                model->add_rbm(std::string("EALL"));
            },
            "InertiaRelief: adjusting external load matrix and adding RBM");
    }

    // (2b) Optional: rebalance loads (independent of inertia relief)
    if (rebalance_loads) {
        logging::error(supps.empty(),
                       "Rebalancing Loads: cannot be used with *SUPPORT in this load case. "
                       "Remove all referenced support collectors.");

        Timer::measure([&]() { fem::rebalance_loads(*model->_data, global_load_mat); },
                       "rebalancing of loads");
    }

    // (3) Constraint equations from supports/ties/couplings
    // Keep grouped constraints to later mask reactions to support-constrained DOFs
    auto groups = Timer::measure(
        [&]() { return this->model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // (4) Active stiffness K with topology scaling/orientation
    auto K = Timer::measure(
        [&]() { return this->model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K(rho^p, theta)"
    );

    // (4a) Reduce loads to active RHS vector f
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f"
    );

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

    // (5) Assemble constraints and build the selected transformation
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::Options copt;
            copt.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,   // maps (node,dof) -> global index or -1
                K.rows(),             // total active DOFs n
                copt
            );
        },
        "building constraint transformer"
    );

    // Diagnostics for the constraint transformer
    logging::info(true           , "");
    logging::info(true           , "Constraint summary");
    logging::up();
    logging::info(true           , "m (rows of C)        : ", CT->report().equations);
    logging::info(true           , "n (cols of C)        : ", CT->report().dofs);
    if (CT->rank_known()) {
        logging::info(true, "rank(C)              : ", CT->rank());
    } else {
        logging::info(true, "rank(C)              : not computed");
    }
    logging::info(true           , "method               : ", CT->method_name());
    logging::info(true           , "solver unknowns       : ", CT->unknowns());
    logging::info(true           , "homogeneous          : ", CT->homogeneous() ? "true" : "false");
    logging::info(true           , "feasible             : ", CT->feasible() ? "true" : "false");
    logging::info(!CT->feasible(), "residual ||C u - d|| : ", CT->report().residual_norm);
    logging::down();

    // Remove temporary RBM again (equations already built)
    if (inertia_relief) {
        logging::error(model->_data != nullptr, "InertiaRelief: model data not initialized");
        logging::error(!model->_data->rbms.empty(),
                       "InertiaRelief: expected a temporary RBM to be present, but rbms is empty");
        model->_data->rbms.pop_back();
    }

    // (6) Reduced system A q = b  with  A = T^T K T,  b = T^T (f - K u_p)
    // Note: Explicit assembly is suitable for direct solvers.
    auto A = Timer::measure(
        [&]() { return CT->assemble_system_matrix(K); },
        "assembling constraint system matrix"
    );
    auto b = Timer::measure(
        [&]() { return CT->assemble_system_rhs(K, f); },
        "assembling constraint system RHS"
    );

    // (6a) Sanity check for NaN/Inf in A and b
    {
        bool badA = false;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<Precision>::InnerIterator it(A, k); it; ++it) {
                if (!std::isfinite(it.value())) {
                    badA = true;
                    break;
                }
            }
            if (badA) break;
        }
        logging::error(!badA, "Matrix A contains NaN/Inf entries");
        logging::error(b.allFinite(), "b contains NaN/Inf entries");
    }

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b, direct_matrix_type); },
        "solving reduced system A q = b"
    );

    // (7) Recover full displacement u and support reactions from multipliers
    auto u = Timer::measure(
        [&]() { return CT->recover_displacement(q); },
        "recovering full displacement vector u"
    );
    auto r_support = Timer::measure(
        [&]() { return CT->support_reactions(K, f, q); },
        "computing support reactions via multipliers (C_supp^T lambda)"
    );

    // (8) Expand vectors to node x 6 matrices for writer output
    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form"
    );
    auto global_react_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r_support); },
        "expanding support reactions to matrix form"
    );

    // (9) Post-processing: nodal stress/strain from displacements
    auto [stress, strain] = Timer::measure(
        [&]() { return model->compute_stress_nodal(global_disp_mat, false); },
        "Interpolating stress and strain at nodes"
    );

    auto shear_flow = Timer::measure(
        [&]() { return model->compute_shear_flow(global_disp_mat); },
        "computing shear-flow output"
    );

    // (10) Topology metrics
    //  - compliance    : element-wise u^T K_e u contributions (as provided by the model)
    //  - dens_grad     : derivative of compliance w.r.t. density (basic SIMP)
    //  - volumes       : element volumes
    //  - angle_grad    : derivative of compliance w.r.t. orientation angles
    auto compliance_raw = model->compute_compliance(global_disp_mat);
    model::Field density_grad   = model->_data->create_field_("DENS_GRAD"     , model::FieldDomain::ELEMENT, 1, false);

    for (Index r = 0; r < static_cast<Index>(model->_data->max_elems); ++r) {
        const Precision rho         = (*density)(r, 0);
        density_grad  (r, 0)        = -exponent * compliance_raw(r, 0) / rho;
    }
    auto volumes        = model->compute_volumes();
    const bool has_orientation = (orientation != nullptr);
    model::Field angle_grad;
    if (has_orientation) {
        angle_grad = model->compute_compliance_angle_derivative(global_disp_mat);
    }

    // (11) Write results
    // Mask reactions to supports only (NaN elsewhere)
    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);
    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 && e.node_id < support_mask.rows() &&
                e.dof < support_mask.cols())
            {
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
    for (Index i = 0; i < reaction_masked.rows; ++i) {
        for (Index j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = global_react_mat(i, j);
            }
        }
    }

    writer->add_loadcase(id, io::writer::WriterStepType::Static);
    writer->write_field(global_disp_mat , "DISPLACEMENT", model->_data.get());
    writer->write_field(strain          , "STRAIN", model->_data.get());
    writer->write_field(stress          , "STRESS", model->_data.get());
    writer->write_field(global_load_mat , "EXTERNAL_FORCES", model->_data.get());
    writer->write_field(reaction_masked , "REACTION_FORCES", model->_data.get());
    writer->write_field(compliance_raw  , "COMPLIANCE", model->_data.get());
    writer->write_field(density_grad    , "DENS_GRAD", model->_data.get());
    writer->write_field(volumes         , "VOLUME", model->_data.get());
    writer->write_field(*density        , "DENSITY", model->_data.get());
    if (shear_flow.rows > 0) {
        writer->write_field(shear_flow, "SHEAR_FLOW", model->_data.get());
    }
    if (has_orientation) {
        writer->write_field(angle_grad  , "ORIENTATION_GRAD", model->_data.get());
        writer->write_field(*orientation, "ORIENTATION", model->_data.get());
    }

    // (13) Diagnostics (optional): projected residual checks
    CT->post_check_static(K, f, q);

    model->_data->element_stiffness_scale = nullptr;
    model->_data->material_orientation    = nullptr;
    model->step_end();
}
}} // namespace fem::loadcase
