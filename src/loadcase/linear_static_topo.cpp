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
 *  (7) Recover u and reactions r = K u - f; expand to node x 6.
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

#include "../constraints/constraint_map.h"
#include "../constraints/constraint_set.h"
#include "../constraints/constraint_transformer.h"
#include "../constraints/equation.h"
#include "../core/logging.h"
#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../solve/eigen.h"

#include <limits>

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
LinearStaticTopo::LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model)
    : LinearStatic(id, writer, model),
      density(model->_data->max_elems, 1),
      orientation(model->_data->max_elems, 3)
{
    density.setOnes();
    orientation.setZero();
}

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
    const ElementData stiffness_scalar = density.array().pow(exponent);
    model->_data->create_data(model::ElementDataEntries::TOPO_STIFFNESS, 1);
    model->_data->create_data(model::ElementDataEntries::TOPO_ANGLES   , 3);
    model->_data->get(model::ElementDataEntries::TOPO_STIFFNESS) = stiffness_scalar;
    model->_data->get(model::ElementDataEntries::TOPO_ANGLES)    = orientation;

    // (1) Unconstrained DOF indexing (node x 6 -> active dof id or -1)
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // (2) Constraint equations from supports/ties/couplings
    // Keep grouped constraints to later mask reactions to support-constrained DOFs
    auto groups = Timer::measure(
        [&]() { return this->model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // (3) Global loads (node x 6) and reduction to active RHS vector f
    auto global_load_mat = Timer::measure(
        [&]() { return this->model->build_load_matrix(loads); },
        "constructing load matrix (node x 6)"
    );
    auto f = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active RHS vector f"
    );

    // (4) Active stiffness K with topology scaling/orientation
    auto K = Timer::measure(
        [&]() { return this->model->build_stiffness_matrix(active_dof_idx_mat, stiffness_scalar); },
        "constructing stiffness matrix K(rho^p, theta)"
    );

    // (5) Build constraint transformer (Set -> Builder -> Map), wrapped for timing
    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::BuildOptions copt;
            // Recommended: scale columns of C for robust QR (rank detection).
            copt.set.scale_columns = true;
            // Optional tolerances (keep defaults unless your constraints are ill-conditioned):
            // copt.builder.rank_tol_rel = 1e-12;
            // copt.builder.feas_tol_rel = 1e-10;
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

    // (6) Reduced system A q = b  with  A = T^T K T,  b = T^T (f - K u_p)
    // Note: Explicit assembly is suitable for direct solvers.
    auto A = Timer::measure(
        [&]() { return CT->assemble_A(K); },
        "assembling A = T^T K T"
    );
    auto b = Timer::measure(
        [&]() { return CT->assemble_b(K, f); },
        "assembling b = T^T (f - K u_p)"
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
        [&]() { return solve(device, method, A, b); },
        "solving reduced system A q = b"
    );

    // (7) Recover full displacement u and reactions r = K u - f
    auto u = Timer::measure(
        [&]() { return CT->recover_u(q); },
        "recovering full displacement vector u"
    );
    auto r = Timer::measure(
        [&]() { return CT->reactions(K, f, q); },
        "computing reactions r = K u - f"
    );

    // (8) Expand vectors to node x 6 matrices for writer output
    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form"
    );
    auto global_force_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r); },
        "expanding reactions to matrix form"
    );

    // (9) Post-processing: nodal stress/strain from displacements
    NodeData stress, strain;
    std::tie(stress, strain) = Timer::measure(
        [&]() { return model->compute_stress_strain(global_disp_mat); },
        "Interpolating stress and strain at nodes"
    );

    // (10) Topology metrics
    //  - compliance_raw: element-wise u^T K_e u contributions (as provided by the model)
    //  - compliance_adj: SIMP-adjusted compliance (rho^p scaling applied)
    //  - dens_grad     : derivative of compliance w.r.t. density (basic SIMP)
    //  - volumes       : element volumes
    //  - angle_grad    : derivative of compliance w.r.t. orientation angles
    ElementData compliance_raw = model->compute_compliance(global_disp_mat);
    ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);
    ElementData volumes        = model->compute_volumes();
    ElementData angle_grad     = model->compute_compliance_angle_derivative(global_disp_mat);

    // (11) Write results
    // Mask reactions to supports only (NaN elsewhere)
    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);
    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 && e.node_id < support_mask.rows() &&
                e.dof >= 0 && e.dof < support_mask.cols())
            {
                if (active_dof_idx_mat(e.node_id, e.dof) != -1) {
                    support_mask(e.node_id, e.dof) = true;
                }
            }
        }
    }
    NodeData reaction_masked = NodeData::Constant(global_force_mat.rows(), global_force_mat.cols(),
                                                  std::numeric_limits<Precision>::quiet_NaN());
    for (int i = 0; i < reaction_masked.rows(); ++i) {
        for (int j = 0; j < reaction_masked.cols(); ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = global_force_mat(i, j);
            }
        }
    }

    writer->add_loadcase(id);
    writer->write_eigen_matrix(global_disp_mat , "DISPLACEMENT");
    writer->write_eigen_matrix(strain          , "STRAIN");
    writer->write_eigen_matrix(stress          , "STRESS");
    writer->write_eigen_matrix(global_load_mat , "EXTERNAL_FORCES");
    writer->write_eigen_matrix(reaction_masked , "REACTION_FORCES");
    writer->write_eigen_matrix(compliance_raw  , "COMPLIANCE_RAW");
    writer->write_eigen_matrix(compliance_adj  , "COMPLIANCE_ADJ");
    writer->write_eigen_matrix(dens_grad       , "DENS_GRAD");
    writer->write_eigen_matrix(volumes         , "VOLUME");
    writer->write_eigen_matrix(density         , "DENSITY");
    writer->write_eigen_matrix(angle_grad      , "ORIENTATION_GRAD");
    writer->write_eigen_matrix(orientation     , "ORIENTATION");

    // (12) Cleanup element scratch data from the model store
    model->_data->remove(model::ElementDataEntries::TOPO_STIFFNESS);
    model->_data->remove(model::ElementDataEntries::TOPO_ANGLES);

    // (13) Diagnostics (optional): projected residual checks
    CT->post_check_static(K, f, u);
}

}} // namespace fem::loadcase
