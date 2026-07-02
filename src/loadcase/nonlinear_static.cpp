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
#include "../writer/write_mtx.h"
#include "tools/regularise_stiffness.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

namespace {

model::Field current_positions_from_displacement(
    const model::Field& reference,
    const model::Field& displacement
) {
    logging::error(reference.domain      == model::FieldDomain::NODE,
                   "NonlinearStatic: reference positions must use NODE domain");
    logging::error(displacement.domain   == model::FieldDomain::NODE,
                   "NonlinearStatic: displacement must use NODE domain");
    logging::error(reference.rows        == displacement.rows,
                   "NonlinearStatic: reference/displacement row mismatch");
    logging::error(reference.components  == displacement.components,
                   "NonlinearStatic: reference/displacement component mismatch");
    logging::error(reference.components  >= 6,
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
    const ConstraintTransformer& transformer,
    const DynamicVector&         external_force,
    const DynamicVector&         internal_force,
    DynamicVector&               reduced_residual
) {
    const DynamicVector residual = external_force - internal_force;
    transformer.apply_Tt(residual, reduced_residual);

    DynamicVector reduced_external;
    transformer.apply_Tt(external_force, reduced_external);

    DynamicVector reduced_internal;
    transformer.apply_Tt(internal_force, reduced_internal);

    const Index reduced_dofs = std::max<Index>(
        static_cast<Index>(reduced_residual.size()), Index(1)
    );

    const Precision inv_sqrt_reduced_dofs =
        Precision(1) / std::sqrt(static_cast<Precision>(reduced_dofs));

    const Precision residual_rms = reduced_residual.norm() * inv_sqrt_reduced_dofs;
    const Precision external_rms = reduced_external.norm() * inv_sqrt_reduced_dofs;
    const Precision internal_rms = reduced_internal.norm() * inv_sqrt_reduced_dofs;

    const Precision denominator = std::max({
        external_rms,
        internal_rms,
        Precision(1)
    });

    return residual_rms / denominator;
}

Precision calculate_arc_length_constraint(
    const DynamicVector& delta_q,
    Precision            delta_lambda,
    Precision            psi2,
    Precision            load_scale2,
    Precision            arc_radius
) {
    return delta_q.dot(delta_q)
         + psi2 * load_scale2 * delta_lambda * delta_lambda
         - arc_radius * arc_radius;
}

Precision calculate_relative_arc_length_residual(
    Precision constraint,
    Precision arc_radius
) {
    return std::abs(constraint) / std::max(
        arc_radius * arc_radius,
        std::numeric_limits<Precision>::epsilon()
    );
}

void build_arc_length_augmented_system(
    const SparseMatrix&  A,
    const DynamicVector& reduced_residual,
    const DynamicVector& reduced_total_load,
    const DynamicVector& delta_q,
    Precision            delta_lambda,
    Precision            constraint,
    Precision            psi2,
    Precision            load_scale2,
    SparseMatrix&        augmented,
    DynamicVector&       augmented_rhs,
    Precision&           load_scale
) {
    const Index n = static_cast<Index>(A.rows());

    load_scale = std::sqrt(load_scale2);

    logging::error(load_scale > Precision(0) && std::isfinite(load_scale),
                   "ArcLength: invalid load scaling in augmented corrector");

    TripletList augmented_triplets;
    augmented_triplets.reserve(static_cast<std::size_t>(A.nonZeros())
                             + static_cast<std::size_t>(2 * n + 1));

    for (Index outer = 0; outer < A.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator entry(A, outer); entry; ++entry) {
            augmented_triplets.emplace_back(entry.row(), entry.col(), entry.value());
        }
    }

    // Use dmu = load_scale * dlambda as the final unknown. This balances the
    // displacement and load parts of the constraint.
    for (Index i = 0; i < n; ++i) {
        const Precision load_entry = -reduced_total_load(i) / load_scale;

        if (load_entry != Precision(0)) {
            augmented_triplets.emplace_back(i, n, load_entry);
        }
    }

    const DynamicVector constraint_q  = Precision(2) * delta_q;
    const Precision     constraint_mu = Precision(2) * psi2 * load_scale * delta_lambda;

    const Precision constraint_scale = std::max({
        constraint_q.template lpNorm<Eigen::Infinity>(),
        std::abs(constraint_mu),
        std::numeric_limits<Precision>::epsilon()
    });

    for (Index i = 0; i < n; ++i) {
        const Precision value = constraint_q(i) / constraint_scale;

        if (value != Precision(0)) {
            augmented_triplets.emplace_back(n, i, value);
        }
    }

    augmented_triplets.emplace_back(n, n, constraint_mu / constraint_scale);

    augmented.resize(n + 1, n + 1);
    augmented.setFromTriplets(augmented_triplets.begin(), augmented_triplets.end());
    augmented.makeCompressed();

    augmented_rhs.resize(n + 1);
    augmented_rhs.head(n) = reduced_residual;
    augmented_rhs(n)      = -constraint / constraint_scale;
}

} // namespace

NonlinearStatic::NonlinearStatic(ID id, reader::Writer* writer, model::Model* model)
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
    logging::error(initial_increment >= minimum_increment && initial_increment <= maximum_increment,
                   "NONLINEARSTATIC requires INITIAL_INCREMENT between MINIMUM_INCREMENT and MAXIMUM_INCREMENT");
    logging::error(max_iterations > 0,
                   "NONLINEARSTATIC requires MAXITER > 0");
    logging::error(tolerance > Precision(0),
                   "NONLINEARSTATIC requires TOL > 0");
    logging::error(arc_length_psi >= Precision(0),
                   "NONLINEARSTATIC requires ARC_LENGTH_PSI >= 0");
    logging::error(arc_length_final_load_threshold >= Precision(0),
                   "NONLINEARSTATIC requires FINAL_LOAD_CONTROL_THRESHOLD >= 0");
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
            ConstraintTransformer::BuildOptions options;
            options.method            = constraint_method;
            options.set.scale_columns = true;

            return std::make_unique<ConstraintTransformer>(
                equations, active_dof_idx_mat, n_active, options
            );
        },
        "building constraint transformer"
    );
    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true                    , "m (rows of C)     : "   , transformer->report().m);
    logging::info(true                    , "n (cols of C)     : "   , transformer->report().n);
    logging::info(true                    , "rank(C)           : "   , transformer->rank());
    logging::info(true                    , "method            : "   , transformer->method_name());
    logging::info(true                    , "masters (n-r)     : "   , transformer->n_master());
    logging::info(true                    , "homogeneous       : "   , transformer->homogeneous() ? "true" : "false");
    logging::info(true                    , "feasible          : "   , transformer->feasible()    ? "true" : "false");
    logging::info(!transformer->feasible(), "residual ||C u - d|| : ", transformer->report().residual_norm);
    logging::down();

    DynamicVector q_total = DynamicVector::Zero(transformer->n_master());
    DynamicVector u_total = transformer->recover_u(q_total);

    DynamicVector reduced_total_load;
    transformer->apply_Tt(f_total, reduced_total_load);

    model::Field displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);

    auto update_positions = [&]() {
        displacement = mattools::expand_vec_to_mat(active_dof_idx_mat, u_total);
        *model->_data->positions = current_positions_from_displacement(reference_positions, displacement);
    };

    update_positions();

    // create fields used below
    model::NodeData internal_mat      {"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    model::NodeData predictor_internal{"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    model::NodeData final_internal    {"INTERNAL_FORCES", model::FieldDomain::NODE, max_nodes, 6};
    internal_mat      .set_zero();
    predictor_internal.set_zero();
    final_internal    .set_zero();

    // settings for the convergence
    const Precision load_increment_growth  = Precision(1.5);
    const Precision load_increment_cutback = Precision(0.5);

    const int fast_convergence_iterations = 4;
    const int slow_convergence_iterations = 10;
    const int maximum_cutbacks            = 20;

    Precision load_factor               = Precision(0);
    Precision load_increment            = initial_increment;
    int       accepted_increment        = 0;
    int       cutback_count             = 0;
    bool      final_load_control_failed = false;

    DynamicVector previous_delta_q      = DynamicVector::Zero(transformer->n_master());
    Precision     previous_delta_lambda = Precision(1);

    logging::info(true, "");
    logging::info(true, "Solver: ", solver::get_solver_name(device, method));
    logging::info(true, "Control: ", control == NonlinearControl::ArcLength ? "ARC LENGTH" : "LOAD CONTROL");
    logging::info(true, "");
    logging::info(true, " inc iter      lambda        rel_res          du_norm   asm_ms solve_ms");
    logging::info(true, "--------------------------------------------------------------------------");

    // Stream converged increment frames as they become available. This keeps
    // the accepted path in the result file even if a later increment fails.
    writer->add_loadcase(id);

    while (load_factor < Precision(1) - tolerance && accepted_increment < max_increments) {
        const bool      arc_length_control      = control == NonlinearControl::ArcLength;
        const Precision remaining_delta_lambda  = Precision(1) - load_factor;
        const bool      final_load_step         =
            arc_length_control &&
            !final_load_control_failed &&
            remaining_delta_lambda > Precision(0) &&
            (remaining_delta_lambda <= arc_length_final_load_threshold ||
             load_increment >= remaining_delta_lambda);
        const bool      arc_length              = arc_length_control && !final_load_step;
        const Precision psi2                    = arc_length_psi * arc_length_psi;

        Precision target_load_factor = load_factor;
        Precision arc_radius         = Precision(0);
        Precision load_scale2        = Precision(0);

        const DynamicVector q_accepted      = q_total;
        const DynamicVector u_accepted      = u_total;
        const Precision     lambda_accepted = load_factor;

        if (!arc_length) {
            if (final_load_step) {
                load_increment     = std::min(load_increment, remaining_delta_lambda);
                target_load_factor = Precision(1);

                logging::info(true,
                              "Trying final load-controlled step from lambda = ", lambda_accepted,
                              " to 1");
            } else {
                load_increment     = std::min(load_increment, remaining_delta_lambda);
                target_load_factor = load_factor + load_increment;
            }
        } else {
            update_positions();

            predictor_internal.set_zero();

            auto Kt_predictor = model->build_tangent_stiffness_matrix(
                active_dof_idx_mat, predictor_internal, displacement
            );

            if (regularize_zero_stiffness_rows) {
                regularise_stiffness(Kt_predictor, zero_stiffness_regularization_alpha);
            }

            auto A_predictor = transformer->assemble_A(Kt_predictor);

            logging::disable();
            DynamicVector dq_load = solver::solve(device, method, A_predictor, reduced_total_load);
            logging::enable();

            load_scale2 = dq_load.squaredNorm();

            logging::error(load_scale2 > Precision(0),
                           "ArcLength: zero load predictor");

            const Precision predictor_norm = std::sqrt(load_scale2 + psi2 * load_scale2);

            Precision path_sign = Precision(1);

            if (accepted_increment > 0) {
                const Precision direction =
                    previous_delta_q.dot(dq_load)
                  + psi2 * load_scale2 * previous_delta_lambda;

                path_sign = direction >= Precision(0) ? Precision(1) : Precision(-1);
            }

            arc_radius = load_increment * predictor_norm;

            q_total     = q_accepted      + path_sign * load_increment * dq_load;
            load_factor = lambda_accepted + path_sign * load_increment;
            u_total     = transformer->recover_u(q_total);

            target_load_factor = load_factor;
        }

        logging::error(load_increment >= minimum_increment,
                       "NONLINEARSTATIC required increment ", load_increment,
                       " is smaller than minimum increment ", minimum_increment);

        bool converged       = false;
        int  iterations_used = 0;

        Timer asm_timer;
        Timer solve_timer;

        // do the iteration for the current load target
        for (int iter = 1; iter <= max_iterations; ++iter) {
            asm_timer.start();

            update_positions();

            // store the internal forces at the current displacement inside a matrix.
            internal_mat.set_zero();

            auto Kt = model->build_tangent_stiffness_matrix(
                active_dof_idx_mat, internal_mat, displacement
            );

            if (regularize_zero_stiffness_rows) {
                regularise_stiffness(Kt, zero_stiffness_regularization_alpha);
            }

            auto f_int = mattools::reduce_mat_to_vec(active_dof_idx_mat, internal_mat);

            const DynamicVector f_ext =
                (arc_length ? load_factor : target_load_factor) * f_total;

            DynamicVector reduced_residual;
            const Precision rel_residual = calculate_relative_residual(
                *transformer, f_ext, f_int, reduced_residual
            );

            Precision rel_arc = Precision(0);

            if (arc_length) {
                const DynamicVector delta_q      = q_total     - q_accepted;
                const Precision     delta_lambda = load_factor - lambda_accepted;

                const Precision arc_constraint = calculate_arc_length_constraint(
                    delta_q, delta_lambda, psi2, load_scale2, arc_radius
                );

                rel_arc = calculate_relative_arc_length_residual(arc_constraint, arc_radius);
            }

            Precision du_norm  = Precision(0);
            Time      solve_ms = Time(0);

            final_internal = internal_mat;

            if (rel_residual <= tolerance && (!arc_length || rel_arc <= tolerance)) {
                asm_timer.stop();

                logging::info(
                    true,
                    std::setw(4), accepted_increment + 1,
                    std::setw(5), iter,
                    std::scientific, std::setprecision(3),
                    std::setw(12), arc_length ? load_factor : target_load_factor,
                    std::setw(15), rel_residual,
                    std::setw(15), du_norm,
                    std::fixed, std::setprecision(1),
                    std::setw(9), asm_timer.elapsed(),
                    std::setw(9), solve_ms
                );

                converged       = true;
                iterations_used = iter;

                break;
            }

            auto A = transformer->assemble_A(Kt);

            logging::error(reduced_residual.allFinite(),
                           "Reduced residual contains NaN/Inf entries");

            asm_timer.stop();

            solve_timer.start();

            DynamicVector dq;
            logging::disable();

            if (!arc_length) {
                dq = solver::solve(device, method, A, reduced_residual);
            } else {
                const DynamicVector delta_q      = q_total     - q_accepted;
                const Precision     delta_lambda = load_factor - lambda_accepted;

                const Precision arc_constraint = calculate_arc_length_constraint(
                    delta_q, delta_lambda, psi2, load_scale2, arc_radius
                );

                // Solve equilibrium and the linearized arc-length constraint
                // as one bordered system. The former partitioned corrector
                // solved K*dq_r=r and K*dq_f=f separately and then selected a
                // root of a quadratic equation. That formulation becomes
                // numerically singular at precisely the load limit points for
                // which arc-length control is needed, even though the bordered
                // continuation system remains regular.
                const Index n = static_cast<Index>(A.rows());

                SparseMatrix  augmented;
                DynamicVector augmented_rhs;
                Precision     load_scale = Precision(0);

                build_arc_length_augmented_system(
                    A, reduced_residual, reduced_total_load,
                    delta_q, delta_lambda, arc_constraint, psi2, load_scale2,
                    augmented, augmented_rhs, load_scale
                );

                const DynamicVector correction = solver::solve(
                    device, method, augmented, augmented_rhs, solver::DirectSolverMatrixType::General
                );

                logging::error(correction.allFinite(),
                               "ArcLength: augmented corrector returned NaN/Inf entries");

                dq           = correction.head(n);
                load_factor += correction(n) / load_scale;
            }

            logging::enable();
            solve_timer.stop();
            solve_ms = solve_timer.elapsed();

            q_total += dq;
            u_total  = transformer->recover_u(q_total);

            DynamicVector du = transformer->map().T() * dq;
            du_norm          = du.norm();

            logging::info(
                true,
                std::setw(4), accepted_increment + 1,
                std::setw(5), iter,
                std::scientific, std::setprecision(3),
                std::setw(12), arc_length ? load_factor : target_load_factor,
                std::setw(15), rel_residual,
                std::setw(15), du_norm,
                std::fixed, std::setprecision(1),
                std::setw(9), asm_timer.elapsed(),
                std::setw(9), solve_ms
            );
        }

        if (!converged) {
            q_total     = q_accepted;
            u_total     = u_accepted;
            load_factor = lambda_accepted;

            update_positions();

            logging::error(adaptive_increments,
                           "NONLINEARSTATIC fixed increment failed at lambda = ",
                           target_load_factor,
                           "; ADAPTIVE=OFF forbids cutback");

            load_increment *= load_increment_cutback;
            cutback_count++;

            logging::info(true,
                          "Increment rejected at lambda = ", target_load_factor,
                          "; reducing increment to ", load_increment);

            logging::error(load_increment >= minimum_increment,
                           "NONLINEARSTATIC required cutback increment is smaller than minimum increment");

            logging::error(cutback_count <= maximum_cutbacks,
                           "NONLINEARSTATIC exceeded maximum number of cutbacks");

            if (final_load_step) {
                final_load_control_failed = true;

                logging::info(true,
                              "Final load-controlled step did not converge; "
                              "continuing with arc-length cutback");
            }

            continue;
        }

        if (!arc_length) {
            load_factor = target_load_factor;
        } else {
            previous_delta_q      = q_total     - q_accepted;
            previous_delta_lambda = load_factor - lambda_accepted;
        }

        accepted_increment          ++;
        cutback_count               = 0;
        final_load_control_failed   = false;

        update_positions();


        auto [increment_stress, increment_strain] =
            model->compute_stress_nodal(displacement, true);
        (void) increment_strain;
        DynamicMatrix lambda_frame(1, 1);
        lambda_frame(0, 0) = load_factor;

        // write displacement and stress for the current increment
        writer->write_field       (displacement    , "DISPLACEMENT_" + std::to_string(accepted_increment), model->_data.get());
        writer->write_field       (increment_stress, "STRESS_"       + std::to_string(accepted_increment), model->_data.get());
        writer->write_eigen_matrix(lambda_frame    , "LAMBDA_"       + std::to_string(accepted_increment));

        if (adaptive_increments) {
            if (iterations_used <= fast_convergence_iterations) {
                load_increment *= load_increment_growth;
            } else if (iterations_used >= slow_convergence_iterations) {
                load_increment *= load_increment_cutback;
            }

            load_increment = std::clamp(load_increment, minimum_increment, maximum_increment);
        } else {
            load_increment = initial_increment;
        }

        logging::info(true,
                      "Accepted increment ", accepted_increment,
                      ": lambda = ", load_factor,
                      ", Newton iterations = ", iterations_used,
                      ", next increment = ", load_increment);
    }

    logging::error(load_factor >= Precision(1) - tolerance,
                   "NONLINEARSTATIC did not reach lambda = 1 within MAX_INCREMENTS");

    update_positions();

    auto [final_stress, final_strain] = Timer::measure(
        [&]() { return model->compute_stress_nodal(displacement, true); },
        "computing final nonlinear nodal stress/strain"
    );

    auto [final_stress_top, final_stress_bot] = Timer::measure(
        [&]() { return model->compute_stress_top_bot(displacement, true); },
        "computing final nonlinear top/bottom stress"
    );

    final_internal.set_zero();

    auto final_Kt = Timer::measure(
        [&]() {
            return model->build_tangent_stiffness_matrix(
                active_dof_idx_mat, final_internal, displacement
            );
        },
        "assembling final nonlinear tangent stiffness K_t and internal force"
    );

    if (!stiffness_file.empty()) {
        if (regularize_zero_stiffness_rows) {
            regularise_stiffness(final_Kt, zero_stiffness_regularization_alpha);
        }

        auto final_A = transformer->assemble_A(final_Kt);

        write_mtx(stiffness_file + "_Kt.mtx", final_Kt);
        write_mtx(stiffness_file + "_A.mtx", final_A);
    }

    auto global_load_final = global_load_total;
    auto reaction_full     = subtract_field(final_internal, global_load_final, "REACTION_FORCES_RAW");

    BooleanMatrix support_mask(active_dof_idx_mat.rows(), active_dof_idx_mat.cols());
    support_mask.setConstant(false);

    for (const auto& eq : groups.supports) {
        for (const auto& e : eq.entries) {
            if (e.node_id >= 0 &&
                e.node_id < support_mask.rows() &&
                e.dof >= 0 &&
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

    for (int i = 0; i < reaction_masked.rows; ++i) {
        for (int j = 0; j < reaction_masked.components; ++j) {
            if (support_mask(i, j)) {
                reaction_masked(i, j) = reaction_full(i, j);
            }
        }
    }

    writer->write_field(displacement     , "DISPLACEMENT"     , model->_data.get());
    writer->write_field(final_strain     , "STRAIN"           , model->_data.get());
    writer->write_field(final_stress     , "STRESS"           , model->_data.get());
    writer->write_field(final_stress_top , "STRESS_TOP"       , model->_data.get());
    writer->write_field(final_stress_bot , "STRESS_BOT"       , model->_data.get());
    writer->write_field(global_load_final, "EXTERNAL_FORCES"  , model->_data.get());
    writer->write_field(final_internal   , "INTERNAL_FORCES"  , model->_data.get());
    writer->write_field(reaction_masked  , "REACTION_FORCES"  , model->_data.get());

    *model->_data->positions = original_positions;
    model->step_end();
}

} // namespace loadcase
} // namespace fem