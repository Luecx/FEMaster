/**
 * @file nonlinear_static.cpp
 * @brief Implements incremental nonlinear static analysis.
 *
 * The nonlinear static load case assembles reduced equilibrium equations for
 * updated nodal positions, applies linear constraints through a null-space
 * transformer, and solves the resulting nonlinear path with either direct load
 * control or arc-length control.
 *
 * Element and contact contributions are evaluated in the current configuration.
 * External loads and prescribed displacements enter through the affine
 * constraint recovery
 *
 *     u(lambda, q) = lambda * u_particular + T q.
 *
 * The load case wires stateful nonlinear subsystems into the generic trial
 * callbacks owned by the path controllers. `NonlinearStateManager` owns the
 * material history buffers and coordinates contact trial state, while the
 * nonlinear solver remains independent of both state implementations.
 *
 * @see NonlinearStatic
 * @see tools::LoadControl
 * @see tools::ArcLengthControl
 * @see tools::NonlinearStateManager
 */

#include "nonlinear_static.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"
#include "../solve/get_solver_name.h"
#include "../io/writer/write_mtx.h"
#include "tools/arc_length_control.h"
#include "tools/load_control.h"
#include "tools/nonlinear_state_manager.h"
#include "tools/regularise_stiffness.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <string>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

/**
 * Builds the current nodal position field from reference coordinates and total
 * generalized displacement.
 *
 * Translational components are interpreted as displacements from the reference
 * configuration and are added to the reference coordinates. Rotational
 * components are generalized displacement coordinates and therefore replace the
 * corresponding entries directly.
 *
 * @param reference Reference nodal coordinates in the undeformed configuration.
 * @param displacement Total generalized nodal displacement field.
 *
 * @return Current nodal position field used by element and contact assembly.
 */
model::Field current_positions_from_displacement(
    const model::Field& reference,
    const model::Field& displacement
) {
    logging::error(reference.domain == model::FieldDomain::NODE,
                   "NonlinearStatic: reference positions must use NODE domain");
    logging::error(displacement.domain == model::FieldDomain::NODE,
                   "NonlinearStatic: displacement must use NODE domain");
    logging::error(reference.rows == displacement.rows,
                   "NonlinearStatic: reference/displacement row mismatch");
    logging::error(reference.components == displacement.components,
                   "NonlinearStatic: reference/displacement component mismatch");
    logging::error(reference.components >= 6,
                   "NonlinearStatic: position fields require six components");

    model::Field current = reference;

    for (Index i = 0; i < current.rows; ++i) {
        for (Index d = 0; d < 3; ++d) {
            current(i, d) = reference(i, d) + displacement(i, d);
        }

        // The rotational entries are the total generalized displacement
        // coordinates rx, ry and rz. They are constrained, solved and written
        // exactly like the translational displacement components.
        for (Index d = 3; d < 6; ++d) {
            current(i, d) = displacement(i, d);
        }
    }

    return current;
}

model::Field subtract_field(
    const model::Field& lhs,
    const model::Field& rhs,
    const std::string&  name
) {
    model::Field result = lhs;
    result.name = name;
    result -= rhs;
    return result;
}

/**
 * Computes the relative residual norm used by nonlinear static convergence.
 *
 * The reduced residual is normalized by the largest RMS scale among the reduced
 * external force, the implied reduced internal force and one. This keeps the
 * convergence measure meaningful for both force-driven and displacement-driven
 * analyses, including states where the external load vector is small.
 *
 * @param reduced_residual Projected nonlinear residual in reduced coordinates.
 * @param reduced_external Projected external load vector at the current load
 *                         factor.
 *
 * @return Dimensionless residual norm used by Newton convergence checks.
 */
Precision calculate_relative_residual(
    const DynamicVector& reduced_residual,
    const DynamicVector& reduced_external
) {
    const DynamicVector reduced_internal =
        reduced_external - reduced_residual;

    const Index reduced_dofs = std::max<Index>(
        static_cast<Index>(reduced_residual.size()),
        Index(1)
    );

    const Precision inv_sqrt_reduced_dofs =
        Precision(1) / std::sqrt(static_cast<Precision>(reduced_dofs));

    const Precision residual_rms =
        reduced_residual.norm() * inv_sqrt_reduced_dofs;
    const Precision external_rms =
        reduced_external.norm() * inv_sqrt_reduced_dofs;
    const Precision internal_rms =
        reduced_internal.norm() * inv_sqrt_reduced_dofs;

    const Precision denominator = std::max({
        external_rms,
        internal_rms,
        Precision(1)
    });

    return residual_rms / denominator;
}

} // namespace

NonlinearStatic::NonlinearStatic(ID id, io::writer::ResultWriters* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

/**
 * Executes the nonlinear static load case.
 *
 * The analysis starts from the reference coordinates and incrementally solves
 * the constrained nonlinear equilibrium equation in reduced coordinates. The
 * constraint transformer represents all active linear constraints as
 *
 *     u = lambda u_p + T q,
 *
 * where `u_p` is the particular solution for prescribed displacements and `T`
 * spans the admissible incremental displacement space.
 *
 * For every residual or tangent evaluation, the supplied reduced state is
 * expanded to total nodal displacement, the current nodal position field is
 * updated, structural element forces are assembled, and nonlinear contact forces
 * are added by the model assembly routine. The residual is then projected back
 * into the reduced space.
 *
 * Stateful nonlinear subsystems are connected through generic trial callbacks:
 *
 * - increment trials cover a complete attempted nonlinear increment and are
 *   committed only after the path controller accepts the increment,
 * - line-search trials cover temporary residual evaluations inside Newton and
 *   are committed only when the step length is accepted,
 * - active-set updates allow discontinuous nonlinear state to be refreshed after
 *   a converged Newton solve and can request a Newton restart at the same load
 *   factor.
 *
 * Solver controls such as maximum increments, maximum iterations and adaptive
 * iteration thresholds are forwarded exactly as configured by the load case. The
 * presence of contact or another stateful subsystem does not change these
 * parameters internally.
 */
void NonlinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "NONLINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Validate path-control, solver and model settings before changing model state
    logging::error(max_increments > 0,
        "NONLINEARSTATIC requires MAX_INCREMENTS > 0");
    logging::error(initial_increment > Precision(0),
        "NONLINEARSTATIC requires INITIAL_INCREMENT > 0");
    logging::error(minimum_increment > Precision(0),
        "NONLINEARSTATIC requires MINIMUM_INCREMENT > 0");
    logging::error(maximum_increment >= minimum_increment,
        "NONLINEARSTATIC requires MAXIMUM_INCREMENT >= MINIMUM_INCREMENT");
    logging::error(initial_increment >= minimum_increment &&
                   initial_increment <= maximum_increment,
        "NONLINEARSTATIC requires INITIAL_INCREMENT between MINIMUM_INCREMENT and MAXIMUM_INCREMENT");
    logging::error(max_iterations > 0,
        "NONLINEARSTATIC requires MAXITER > 0");
    logging::error(tolerance > Precision(0),
        "NONLINEARSTATIC requires TOL > 0");
    logging::error(arc_length_psi >= Precision(0),
        "NONLINEARSTATIC requires ARC_LENGTH_PSI >= 0");
    logging::error(growth_factor > Precision(0),
        "NONLINEARSTATIC requires GROWTH_FACTOR > 0");
    logging::error(cutback_factor > Precision(0) && cutback_factor < Precision(1),
        "NONLINEARSTATIC requires CUTBACK_FACTOR between 0 and 1");
    logging::error(fast_iterations > 0,
        "NONLINEARSTATIC requires FAST_ITERATIONS > 0");
    logging::error(slow_iterations >= fast_iterations,
        "NONLINEARSTATIC requires SLOW_ITERATIONS >= FAST_ITERATIONS");
    logging::error(maximum_cutbacks > 0,
        "NONLINEARSTATIC requires MAXIMUM_CUTBACKS > 0");
    logging::error(constraint_method == ConstraintTransformer::Method::NullSpace,
        "NONLINEARSTATIC currently supports only NULLSPACE constraints");
    logging::error(method == solver::DIRECT,
        "NONLINEARSTATIC currently supports only DIRECT solver method");
    logging::error(model->_data->positions != nullptr,
        "NonlinearStatic: positions field not initialized");
    logging::error(model->_data->positions_reference != nullptr,
        "NonlinearStatic: positions_reference field not initialized");

    // Preserve the caller-visible position field and start the nonlinear step
    // from the undeformed reference configuration
    const model::Field original_positions  = *model->_data->positions;
    const model::Field reference_positions = *model->_data->positions_reference;
    *model->_data->positions = reference_positions;

    // Bind element sections and initialize step-local element geometry/state caches
    model->assign_sections();
    model->step_begin();

    // Own material history buffers and coordinate all contact trial states for
    // the duration of this nonlinear solution.
    tools::NonlinearStateManager nonlinear_state(*model);

    // Build active DOFs, external loading and constraint groups shared by all
    // nonlinear evaluations
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    auto global_load_total = Timer::measure(
        [&]() { return model->build_load_matrix(loads); },
        "constructing total load matrix (node x 6)"
    );

    auto f_total = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_total); },
        "reducing total load matrix -> active RHS vector f"
    );

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat, supps); },
        "building constraints"
    );

    // Report the semantic constraint groups before flattening them for algebraic use
    report_constraint_groups(groups);

    // Flatten semantic groups into the complete equation set used by the transformer
    auto equations = groups.flatten();

    // Cache global dimensions required by reduced-vector and nodal-field assembly
    const Index n_active  = active_dof_idx_mat.maxCoeff() + 1;
    const Index max_nodes = static_cast<Index>(model->_data->max_nodes);

    // Construct the affine null-space transformation for prescribed and coupled DOFs
    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;

            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                n_active,
                options
            );
        },
        "building constraint transformer"
    );

    logging::info(true                    , "");
    logging::info(true                    , "Constraint summary");
    logging::up();
    logging::info(true                    , "m (rows of C)     : "   , transformer->report().equations);
    logging::info(true                    , "n (cols of C)     : "   , transformer->report().dofs);
    logging::info(true                    , "rank(C)           : "   , transformer->rank());
    logging::info(true                    , "method            : "   , transformer->method_name());
    logging::info(true                    , "solver unknowns   : "   , transformer->unknowns());
    logging::info(true                    , "homogeneous       : "   , transformer->homogeneous() ? "true" : "false");
    logging::info(true                    , "feasible          : "   , transformer->feasible()    ? "true" : "false");
    logging::info(!transformer->feasible(), "residual ||C u - d|| : ", transformer->report().residual_norm);
    logging::down();

    DynamicVector       q_total      = DynamicVector::Zero(transformer->unknowns());
    const DynamicVector u_particular = transformer->recover_displacement(q_total);

    logging::error(control != NonlinearControl::ArcLength || transformer->homogeneous(),
        "NONLINEARSTATIC ARC LENGTH does not support prescribed displacements");

    auto recover_total_displacement = [&](const DynamicVector& q, Precision lambda) -> DynamicVector {
        return lambda * u_particular + transformer->recover_increment(q);
    };

    DynamicVector u_total = recover_total_displacement(q_total, Precision(0));

    DynamicVector reduced_total_load;
    transformer->project_vector(f_total, reduced_total_load);

    model::Field displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);

    auto update_positions = [&]() {
        displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);
        *model->_data->positions = current_positions_from_displacement(
            reference_positions,
            displacement
        );
    };

    update_positions();

    model::Field final_internal{"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    final_internal.set_zero();

    Precision load_factor = Precision(0);

    logging::info(true, "");
    logging::info(true, "Solver: ", solver::get_solver_name(device, method));
    logging::info(true, "Control: ",
        control == NonlinearControl::ArcLength ? "ARC LENGTH" : "LOAD CONTROL");
    logging::info(true, "");
    logging::info(true, " inc iter      lambda        rel_res          du_norm   ls   asm_ms solve_ms");
    logging::info(true, "----------------------------------------------------------------------------");

    // Stream converged increment frames as they become available. This keeps
    // the accepted path in the result file even if a later increment fails.
    writer->add_loadcase(id, io::writer::WriterStepType::Static);

    Index last_converged_increment = 0;

    auto assemble_state = [&](const DynamicVector& q,
                              Precision            lambda,
                              DynamicVector&       residual,
                              SparseMatrix&        tangent,
                              SparseMatrix*        full_tangent) {
        // Every independent nonlinear evaluation starts from the material
        // history of the last accepted physical increment.
        nonlinear_state.reset_material_state();

        // Evaluate the supplied solver state without modifying q_total. This is
        // essential because line-search trial states must not modify the actual
        // Newton iterate.
        const DynamicVector u_evaluation =
            recover_total_displacement(q, lambda);

        const model::Field displacement_evaluation =
            mattools::expand_vec_to_mat(active_dof_idx_mat, u_evaluation);

        *model->_data->positions = current_positions_from_displacement(
            reference_positions,
            displacement_evaluation
        );

        model::NodeData internal_mat{
            "INTERNAL_FORCES",
            model::FieldDomain::NODE,
            max_nodes,
            6
        };
        internal_mat.set_zero();

        SparseMatrix Kt;

        const bool logging_was_enabled = logging::is_enabled();
        logging::disable();

        try {
            Kt = model->build_tangent_stiffness_matrix(
                active_dof_idx_mat,
                internal_mat,
                displacement_evaluation
            );
        } catch (...) {
            if (logging_was_enabled) {
                logging::enable();
            }

            throw;
        }

        if (logging_was_enabled) {
            logging::enable();
        }

        if (regularize_zero_stiffness_rows) {
            regularise_stiffness(
                Kt,
                zero_stiffness_regularization_alpha
            );
        }

        const DynamicVector internal_force =
            mattools::reduce_mat_to_vec(
                active_dof_idx_mat,
                internal_mat
            );

        const DynamicVector external_force = lambda * f_total;
        const DynamicVector full_residual  = external_force - internal_force;

        transformer->project_vector(full_residual, residual);

        if (full_tangent) {
            *full_tangent = Kt;
        }

        tangent = transformer->assemble_system_matrix(Kt);

        final_internal = internal_mat;

        logging::error(
            residual.allFinite(),
            "Reduced residual contains NaN/Inf entries"
        );
    };

    auto evaluate = [&](const DynamicVector& q,
                        Precision            lambda,
                        DynamicVector&       residual,
                        SparseMatrix&        tangent) {
        assemble_state(
            q,
            lambda,
            residual,
            tangent,
            nullptr
        );
    };

    auto evaluate_residual = [&](const DynamicVector& q,
                                 Precision            lambda,
                                 DynamicVector&       residual) {
        // Line-search candidates are independent material evaluations of the
        // same physical increment and therefore restart from committed history.
        nonlinear_state.reset_material_state();

        // Evaluate the same nonlinear residual as the full tangent assembly,
        // but avoid constructing the material and geometric element tangents.
        // Line-search acceptance only uses the projected residual norm.
        const DynamicVector u_evaluation =
            recover_total_displacement(q, lambda);

        const model::Field displacement_evaluation =
            mattools::expand_vec_to_mat(active_dof_idx_mat, u_evaluation);

        *model->_data->positions = current_positions_from_displacement(
            reference_positions,
            displacement_evaluation
        );

        model::NodeData internal_mat{
            "INTERNAL_FORCES",
            model::FieldDomain::NODE,
            max_nodes,
            6
        };
        internal_mat.set_zero();

        const bool logging_was_enabled = logging::is_enabled();
        logging::disable();

        try {
            model->build_internal_force_nonlinear(
                active_dof_idx_mat,
                internal_mat,
                displacement_evaluation
            );
        } catch (...) {
            if (logging_was_enabled) {
                logging::enable();
            }

            throw;
        }

        if (logging_was_enabled) {
            logging::enable();
        }

        const DynamicVector internal_force =
            mattools::reduce_mat_to_vec(
                active_dof_idx_mat,
                internal_mat
            );

        const DynamicVector external_force = lambda * f_total;
        const DynamicVector full_residual  = external_force - internal_force;

        transformer->project_vector(full_residual, residual);

        final_internal = internal_mat;

        logging::error(residual.allFinite(),
            "Reduced residual contains NaN/Inf entries");
    };

    auto linear_solve = [&](const SparseMatrix&  tangent,
                            const DynamicVector& rhs) {
        SparseMatrix matrix = tangent;

        DynamicVector solution;
        const bool logging_was_enabled = logging::is_enabled();
        logging::disable();

        try {
            solution = solver::solve(device, method, matrix, rhs, solver::DirectSolverMatrixType::General);
        } catch (...) {
            if (logging_was_enabled) {
                logging::enable();
            }
            throw;
        }

        if (logging_was_enabled) {
            logging::enable();
        }
        return solution;
    };

    auto predictor = [&](DynamicVector& q,
                         Precision      lambda,
                         Precision      target_lambda) {
        if (q.size() == 0) {
            return;
        }

        /*
         * Assemble the tangent at the last accepted equilibrium state. For
         *
         *     u(lambda, q) = lambda * u_p + T q
         *
         * the equilibrium-path tangent follows from
         *
         *     T^T K_t T dq/dlambda = T^T (f - K_t u_p).
         *
         * ConstraintTransformer::assemble_system_rhs() already evaluates the
         * right-hand side T^T (f - K_t u_p), so the same predictor handles
         * force loading, prescribed displacements and arbitrary mixtures.
         */
        DynamicVector accepted_residual;
        SparseMatrix  accepted_tangent;
        SparseMatrix  accepted_full_tangent;

        /*
         * The outer increment trial is deliberately still unfrozen here. The
         * predictor tangent is evaluated at the last accepted configuration,
         * but it must not select and freeze the partners that Newton will use
         * at the predicted target configuration. A nested frozen trial makes
         * this assembly completely state-neutral and is always rolled back.
         */
        nonlinear_state.begin_contact_frozen_trial();

        try {
            assemble_state(
                q,
                lambda,
                accepted_residual,
                accepted_tangent,
                &accepted_full_tangent
            );
        } catch (...) {
            nonlinear_state.rollback_contact_trial();
            throw;
        }

        nonlinear_state.rollback_contact_trial();

        const DynamicVector predictor_rhs =
            transformer->assemble_system_rhs(
                accepted_full_tangent,
                f_total
            );

        const DynamicVector dq_dlambda =
            linear_solve(
                accepted_tangent,
                predictor_rhs
            );

        q += (target_lambda - lambda) * dq_dlambda;
    };

    auto matrix_solve = [&](const SparseMatrix&  tangent,
                            const DynamicMatrix& rhs) {
        SparseMatrix matrix = tangent;

        DynamicMatrix solution;
        const bool logging_was_enabled = logging::is_enabled();
        logging::disable();

        try {
            solution = solver::solve(
                device,
                method,
                matrix,
                rhs,
                solver::DirectSolverMatrixType::General
            );
        } catch (...) {
            if (logging_was_enabled) {
                logging::enable();
            }
            throw;
        }

        if (logging_was_enabled) {
            logging::enable();
        }
        return solution;
    };

    auto residual_norm = [&](const DynamicVector& residual,
                             Precision            lambda) {
        const DynamicVector reduced_external = lambda * reduced_total_load;

        return calculate_relative_residual(residual, reduced_external);
    };

    auto correction_norm = [&](const DynamicVector& q,
                               const DynamicVector& dq) {
        (void) q;

        const DynamicVector du = transformer->recover_increment(dq);
        return du.norm();
    };

    auto on_iteration = [&](Index     increment,
                            Index     iteration,
                            Precision lambda,
                            Precision residual_norm,
                            Precision correction_norm,
                            Precision convergence_order,
                            Index     line_search_iterations,
                            Time      assembly_ms,
                            Time      solve_ms,
                            bool      converged) {
        (void) convergence_order;
        (void) converged;

        logging::info(
            true,
            std::setw(4), increment,
            std::setw(5), iteration,
            std::scientific, std::setprecision(3),
            std::setw(12), lambda,
            std::setw(15), residual_norm,
            std::setw(15), correction_norm,
            std::setw(5), line_search_iterations,
            std::fixed, std::setprecision(1),
            std::setw(9), assembly_ms,
            std::setw(9), solve_ms
        );
    };

    auto on_increment = [&](Index                increment,
                            const DynamicVector& q,
                            Precision            lambda) {
        last_converged_increment = increment;

        q_total     = q;
        load_factor = lambda;
        u_total     = recover_total_displacement(q_total, lambda);

        update_positions();

        writer->write_field(
            displacement,
            "DISPLACEMENT_" + std::to_string(increment),
            model->_data.get(),
            lambda
        );

        auto [increment_stress, increment_strain] =
            model->compute_stress_nodal(displacement, true);
        (void) increment_strain;

        writer->write_field(
            increment_stress,
            "STRESS_" + std::to_string(increment),
            model->_data.get(),
            lambda
        );

        model::Field lambda_field{"LAMBDA_" + std::to_string(increment), model::FieldDomain::UNKNOWN, 1, 1};
        lambda_field(0) = lambda;

        writer->write_field(lambda_field, lambda_field.name, nullptr);
    };

    // Start a complete increment attempt. Material history is reset to the last
    // accepted increment and contact may update partner ownership once before
    // freezing the discrete contact problem for Newton.
    auto begin_increment_trial = [&]() {
        nonlinear_state.reset_material_state();
        nonlinear_state.begin_contact_update_trial();
    };

    // Accept the converged material history and the outer contact trial.
    auto commit_increment_trial = [&]() {
        nonlinear_state.commit_contact_trial();
        nonlinear_state.commit_material_state();
    };

    // Reject the complete increment attempt and restore both state subsystems to
    // the last accepted physical increment.
    auto rollback_increment_trial = [&]() {
        nonlinear_state.rollback_contact_trial();
        nonlinear_state.reset_material_state();
    };

    // Open a temporary contact trial for one line-search candidate. Material
    // history is reset by evaluate_residual() immediately before assembly.
    auto begin_line_search_trial = [&]() {
        nonlinear_state.begin_contact_frozen_trial();
    };

    // Promote the accepted candidate's contact geometry to its parent trial
    auto commit_line_search_trial = [&]() {
        nonlinear_state.commit_contact_trial();
    };

    // Discard a rejected candidate without changing its parent contact state
    auto rollback_line_search_trial = [&]() {
        nonlinear_state.rollback_contact_trial();
    };

    // Refresh discontinuous contact state after Newton convergence. The update
    // trial evaluates the converged geometry with partner updates enabled once;
    // a changed partner signature or multiplier requests another Newton solve.
    auto update_active_set = [&](const DynamicVector& q,
                                 Precision            lambda) {
        if (model->_data->contacts.empty()) {
            return true;
        }

        nonlinear_state.begin_contact_update_trial();

        try {
            DynamicVector active_set_residual;
            SparseMatrix  active_set_tangent;

            assemble_state(
                q,
                lambda,
                active_set_residual,
                active_set_tangent,
                nullptr
            );
        } catch (...) {
            nonlinear_state.rollback_contact_trial();
            throw;
        }

        const bool unchanged = nonlinear_state.update_contact_active_set();
        nonlinear_state.commit_contact_trial();
        return unchanged;
    };

    bool        converged      = false;
    std::string failure_reason = "NONE";

    // Use the nonlinear controls exactly as configured. The trial callbacks
    // below manage stateful nonlinear subsystems without changing increment or
    // iteration limits based on which subsystem is present.
    const Index configured_max_iterations =
        static_cast<Index>(max_iterations);

    const Index configured_max_increments =
        static_cast<Index>(max_increments);

    const Index configured_slow_iterations =
        static_cast<Index>(slow_iterations);

    if (control == NonlinearControl::LoadControl) {
        tools::LoadControl load_control;

        load_control.maximum_increments = configured_max_increments;
        load_control.maximum_iterations = configured_max_iterations;
        load_control.tolerance          = tolerance;
        load_control.initial_increment  = initial_increment;
        load_control.minimum_increment  = minimum_increment;
        load_control.maximum_increment  = maximum_increment;
        load_control.growth_factor      = growth_factor;
        load_control.cutback_factor     = cutback_factor;
        load_control.fast_iterations    = static_cast<Index>(fast_iterations);
        load_control.slow_iterations    = configured_slow_iterations;
        load_control.maximum_cutbacks   = static_cast<Index>(maximum_cutbacks);
        load_control.adaptive           = adaptive_increments;

        load_control.begin_increment_trial    = begin_increment_trial;
        load_control.commit_increment_trial   = commit_increment_trial;
        load_control.rollback_increment_trial = rollback_increment_trial;

        load_control.begin_line_search_trial    = begin_line_search_trial;
        load_control.commit_line_search_trial   = commit_line_search_trial;
        load_control.rollback_line_search_trial = rollback_line_search_trial;

        load_control.update_active_set = update_active_set;

        converged = load_control.solve(
            q_total,
            load_factor,
            evaluate,
            linear_solve,
            residual_norm,
            correction_norm,
            on_iteration,
            on_increment,
            predictor,
            evaluate_residual
        );

        failure_reason = load_control.failure_reason();
    } else {
        tools::ArcLengthControl arc_length_control;

        arc_length_control.maximum_increments = configured_max_increments;
        arc_length_control.maximum_iterations = configured_max_iterations;
        arc_length_control.tolerance          = tolerance;
        arc_length_control.initial_increment  = initial_increment;
        arc_length_control.minimum_increment  = minimum_increment;
        arc_length_control.maximum_increment  = maximum_increment;
        arc_length_control.psi                = arc_length_psi;
        arc_length_control.growth_factor      = growth_factor;
        arc_length_control.cutback_factor     = cutback_factor;
        arc_length_control.fast_iterations    = static_cast<Index>(fast_iterations);
        arc_length_control.slow_iterations    = configured_slow_iterations;
        arc_length_control.maximum_cutbacks   = static_cast<Index>(maximum_cutbacks);
        arc_length_control.adaptive           = adaptive_increments;

        arc_length_control.begin_increment_trial    = begin_increment_trial;
        arc_length_control.commit_increment_trial   = commit_increment_trial;
        arc_length_control.rollback_increment_trial = rollback_increment_trial;

        arc_length_control.begin_line_search_trial    = begin_line_search_trial;
        arc_length_control.commit_line_search_trial   = commit_line_search_trial;
        arc_length_control.rollback_line_search_trial = rollback_line_search_trial;

        arc_length_control.update_active_set = update_active_set;

        converged = arc_length_control.solve(
            q_total,
            load_factor,
            reduced_total_load,
            evaluate,
            linear_solve,
            matrix_solve,
            residual_norm,
            correction_norm,
            on_iteration,
            on_increment,
            evaluate_residual
        );

        failure_reason = arc_length_control.failure_reason();
    }

    logging::error(converged,
        "NONLINEARSTATIC failed: ", failure_reason);

    // Restore the accepted configuration before independent final result evaluations
    update_positions();

    // Each postprocessing path starts from committed material history because
    // result recovery may call an in-place constitutive model
    nonlinear_state.reset_material_state();
    auto [final_stress, final_strain] = Timer::measure(
        [&]() { return model->compute_stress_nodal(displacement, true); },
        "computing final nonlinear nodal stress/strain"
    );

    // Top/bottom shell recovery is independent of the preceding nodal recovery
    nonlinear_state.reset_material_state();
    auto [final_stress_top, final_stress_bot] = Timer::measure(
        [&]() { return model->compute_stress_top_bot(displacement, true); },
        "computing final nonlinear top/bottom stress"
    );

    final_internal = model::NodeData{"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    final_internal.set_zero();

    // Reassemble the final tangent and matching internal force from committed history
    nonlinear_state.reset_material_state();
    auto final_Kt = Timer::measure(
        [&]() {
            return model->build_tangent_stiffness_matrix(
                active_dof_idx_mat,
                final_internal,
                displacement
            );
        },
        "assembling final nonlinear tangent stiffness K_t and internal force"
    );

    if (!stiffness_file.empty()) {
        if (regularize_zero_stiffness_rows) {
            regularise_stiffness(final_Kt, zero_stiffness_regularization_alpha);
        }

        auto final_A = transformer->assemble_system_matrix(final_Kt);

        io::writer::write_mtx(stiffness_file + "_Kt.mtx", final_Kt);
        io::writer::write_mtx(stiffness_file + "_A.mtx", final_A);
    }

    auto global_load_final = global_load_total;
    global_load_final *= load_factor;
    auto reaction_full     = subtract_field(
        final_internal,
        global_load_final,
        "REACTION_FORCES_RAW"
    );

    BooleanMatrix support_mask(
        active_dof_idx_mat.rows(),
        active_dof_idx_mat.cols()
    );
    support_mask.setConstant(false);

    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 &&
                e.node_id < support_mask.rows() &&
                e.dof < support_mask.cols() &&
                active_dof_idx_mat(e.node_id, e.dof) != -1) {
                support_mask(e.node_id, e.dof) = true;
            }
        }
    }

    model::Field reaction_masked{
        "REACTION_FORCES",
        model::FieldDomain::NODE,
        reaction_full.rows,
        reaction_full.components
    };
    reaction_masked.fill_nan();

    for (Index i = 0; i < reaction_masked.rows; ++i) {
        for (Index j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = reaction_full(i, j);
            }
        }
    }

    const Index final_frame = last_converged_increment > 0
        ? last_converged_increment + 1
        : 1;

    const std::string suffix = "_" + std::to_string(final_frame);

    writer->write_field(displacement     , "DISPLACEMENT"    + suffix, model->_data.get(), load_factor);
    writer->write_field(final_strain     , "STRAIN"          + suffix, model->_data.get(), load_factor);
    writer->write_field(final_stress     , "STRESS"          + suffix, model->_data.get(), load_factor);
    writer->write_field(final_stress_top , "STRESS_TOP"      + suffix, model->_data.get(), load_factor);
    writer->write_field(final_stress_bot , "STRESS_BOT"      + suffix, model->_data.get(), load_factor);
    writer->write_field(global_load_final, "EXTERNAL_FORCES" + suffix, model->_data.get(), load_factor);
    writer->write_field(final_internal   , "INTERNAL_FORCES" + suffix, model->_data.get(), load_factor);
    writer->write_field(reaction_masked  , "REACTION_FORCES" + suffix, model->_data.get(), load_factor);

    // Restore the caller-visible model configuration and release step-local caches
    *model->_data->positions = original_positions;
    model->step_end();
}

} // namespace loadcase
} // namespace fem
