/**
 * @file nonlinear_static.cpp
 * @brief Implements incremental nonlinear static analysis.
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
#include "tools/regularise_stiffness.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

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

void NonlinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "NONLINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // perform important checks
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

    // create a copy of original and reference positions and reset them at the end
    const model::Field original_positions  = *model->_data->positions;
    const model::Field reference_positions = *model->_data->positions_reference;
    *model->_data->positions = reference_positions;

    // assign sections so every element knows its section and the elements can prepare
    model->assign_sections();
    model->step_begin();

    // compute general stuff needed below
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

    // report groups if desired by the user
    report_constraint_groups(groups);

    // convert into actual equations (just flattening all the equation groups)
    auto equations = groups.flatten();

    // track a few simple coefficients used below
    const Index n_active  = active_dof_idx_mat.maxCoeff() + 1;
    const Index max_nodes = static_cast<Index>(model->_data->max_nodes);

    // create the constraint transformer
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
    logging::info(true, " inc iter      lambda        rel_res          du_norm   asm_ms solve_ms");
    logging::info(true, "--------------------------------------------------------------------------");

    // Stream converged increment frames as they become available. This keeps
    // the accepted path in the result file even if a later increment fails.
    writer->add_loadcase(id, io::writer::WriterStepType::Static);

    Index last_converged_increment = 0;

    // Initialize the load-control predictor history
    DynamicVector q_previous = q_total;
    DynamicVector q_accepted = q_total;

    Precision lambda_previous = Precision(0);
    Precision lambda_accepted = Precision(0);

    bool secant_available = false;

    auto evaluate = [&](const DynamicVector& q,
                        Precision            lambda,
                        DynamicVector&       residual,
                        SparseMatrix&        tangent) {
        q_total = q;
        u_total = recover_total_displacement(q_total, lambda);
        update_positions();

        model::NodeData internal_mat{
            "INTERNAL_FORCES",
            model::FieldDomain::NODE,
            max_nodes,
            6
        };
        internal_mat.set_zero();

        auto Kt = model->build_tangent_stiffness_matrix(
            active_dof_idx_mat,
            internal_mat,
            displacement
        );

        if (regularize_zero_stiffness_rows) {
            regularise_stiffness(Kt, zero_stiffness_regularization_alpha);
        }

        const DynamicVector internal_force = mattools::reduce_mat_to_vec(
            active_dof_idx_mat,
            internal_mat
        );

        const DynamicVector external_force = lambda * f_total;
        const DynamicVector full_residual  = external_force - internal_force;

        transformer->project_vector(full_residual, residual);
        tangent = transformer->assemble_system_matrix(Kt);
        final_internal = internal_mat;

        logging::error(residual.allFinite(),
            "Reduced residual contains NaN/Inf entries");
    };

    auto linear_solve = [&](const SparseMatrix&  tangent,
                            const DynamicVector& rhs) {
        SparseMatrix matrix = tangent;

        logging::disable();
        DynamicVector solution = solver::solve(device, method, matrix, rhs);
        logging::enable();
        return solution;
    };

    auto predictor = [&](DynamicVector& q,
                         Precision      lambda,
                         Precision      target_lambda) {
        if (!secant_available) {
            return;
        }

        // Extrapolate the last converged reduced increment
        const Precision previous_increment = lambda - lambda_previous;
        const Precision target_increment   = target_lambda - lambda;

        logging::error(previous_increment > Precision(0),
            "NONLINEARSTATIC predictor requires a positive previous increment");

        q = q_accepted + target_increment / previous_increment * (q_accepted - q_previous);
    };

    auto matrix_solve = [&](const SparseMatrix&  tangent,
                            const DynamicMatrix& rhs) {
        SparseMatrix matrix = tangent;

        logging::disable();
        DynamicMatrix solution = solver::solve(
            device,
            method,
            matrix,
            rhs,
            solver::DirectSolverMatrixType::General
        );
        logging::enable();
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
            std::fixed, std::setprecision(1),
            std::setw(9), assembly_ms,
            std::setw(9), solve_ms
        );
    };

    auto on_increment = [&](Index                increment,
                            const DynamicVector& q,
                            Precision            lambda) {
        last_converged_increment = increment;

        q_previous       = q_accepted;
        q_accepted       = q;
        lambda_previous  = lambda_accepted;
        lambda_accepted  = lambda;
        secant_available = true;

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

    bool        converged      = false;
    const char* failure_reason = "NONE";

    if (control == NonlinearControl::LoadControl) {
        tools::LoadControl load_control;

        load_control.maximum_increments = static_cast<Index>(max_increments);
        load_control.maximum_iterations = static_cast<Index>(max_iterations);
        load_control.tolerance          = tolerance;
        load_control.initial_increment  = initial_increment;
        load_control.minimum_increment  = minimum_increment;
        load_control.maximum_increment  = maximum_increment;
        load_control.growth_factor      = growth_factor;
        load_control.cutback_factor     = cutback_factor;
        load_control.fast_iterations    = static_cast<Index>(fast_iterations);
        load_control.slow_iterations    = static_cast<Index>(slow_iterations);
        load_control.maximum_cutbacks   = static_cast<Index>(maximum_cutbacks);
        load_control.adaptive           = adaptive_increments;

        converged = load_control.solve(
            q_total,
            load_factor,
            evaluate,
            linear_solve,
            residual_norm,
            correction_norm,
            on_iteration,
            on_increment,
            predictor
        );

        failure_reason = load_control.failure_reason();
    } else {
        tools::ArcLengthControl arc_length_control;

        arc_length_control.maximum_increments = static_cast<Index>(max_increments);
        arc_length_control.maximum_iterations = static_cast<Index>(max_iterations);
        arc_length_control.tolerance          = tolerance;
        arc_length_control.initial_increment  = initial_increment;
        arc_length_control.minimum_increment  = minimum_increment;
        arc_length_control.maximum_increment  = maximum_increment;
        arc_length_control.psi                = arc_length_psi;
        arc_length_control.growth_factor      = growth_factor;
        arc_length_control.cutback_factor     = cutback_factor;
        arc_length_control.fast_iterations    = static_cast<Index>(fast_iterations);
        arc_length_control.slow_iterations    = static_cast<Index>(slow_iterations);
        arc_length_control.maximum_cutbacks   = static_cast<Index>(maximum_cutbacks);
        arc_length_control.adaptive           = adaptive_increments;

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
            on_increment
        );

        failure_reason = arc_length_control.failure_reason();
    }

    logging::error(converged,
        "NONLINEARSTATIC failed: ", failure_reason);

    update_positions();

    auto [final_stress, final_strain] = Timer::measure(
        [&]() { return model->compute_stress_nodal(displacement, true); },
        "computing final nonlinear nodal stress/strain"
    );

    auto [final_stress_top, final_stress_bot] = Timer::measure(
        [&]() { return model->compute_stress_top_bot(displacement, true); },
        "computing final nonlinear top/bottom stress"
    );

    final_internal = model::NodeData{"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    final_internal.set_zero();

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

    *model->_data->positions = original_positions;
    model->step_end();
}

} // namespace loadcase
} // namespace fem
