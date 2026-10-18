/**
 * @file linear_static.cpp
 * @brief Implements the linear static load case leveraging constraint maps.
 *
 * The structural analysis reads loads and support prescriptions directly from
 * ModelData, reduces K_T du = delta f through the chosen constraint map and recovers
 * global response fields. Optional inertia relief and load rebalancing use the
 * actual active support state. Condition replacement semantics belong to the DSL;
 * model assembly owns physical load integration and equation construction.
 *
 * @see Model::build_load_matrix
 * @see Model::collect_constraints
 * @see constraint::ConstraintTransformer
 */

#include "linear_static.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/rbm.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/all_finite.h"
#include "../mattools/mask_field.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/element/element_structural.h"
#include "../model/model.h"
#include "../solve/eigval/solve_eigval.h"
#include "../io/writer/write_mtx.h"
#include "tools/inertia_relief.h"
#include "tools/rebalance_loads.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>

namespace fem {
namespace loadcase {

using constraint::ConstraintTransformer;

/**
 * Executes constrained linear equilibrium against the current model state.
 *
 * Sections and element caches are prepared before enumerating global DOFs.
 * Direct structural history and collector-activated conditions supply f and C u = d.
 * Optional inertia relief or load rebalancing requires empty structural support
 * state; balancing changes only the assembled force field, and the temporary
 * inertia-relief RBM is removed after constraint collection.
 *
 * Every linear static step solves for an increment around its prepared base
 * configuration. Resume selects the converged base; independent steps use
 * the reference configuration. Only perturbations discard the increment.
 *
 * The affine constraint map reduces K_T du = delta f to the configured backend.
 * Expanded global displacements drive stress/strain and reaction recovery, then
 * requested result fields are written and analysis-local element caches released.
 * Model condition definitions and collector-derived activations are read, never copied into
 * the analysis or modified by this procedure.
 */
void LinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Prepare section assignments and material history before element assembly.
    model->assign_sections();
    model->step_begin();
    if (perturbation) {
        logging::error(model->_data->contacts.empty(),
            "LINEARSTATIC PERTURBATION: contact tangent is not supported");
        logging::error(!inertia_relief && !rebalance_loads,
            "LINEARSTATIC PERTURBATION: inertia relief and rebalancing are not supported");
        logging::error(method == solver::DIRECT,
            "LINEARSTATIC PERTURBATION: a DIRECT solver is required for general tangents");
    }

    // Structural temperatures are prescribed fields, not thermal constraints.
    // The solver evaluates thermal forces at the step's target temperature.
    model::Field* target_temperature = model->_data->temperature_target
        ? model->_data->temperature_target.get()
        : model->_data->temperature.get();

    // Enumerate the current structural system in global generalized coordinates
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_structural_dof_index_matrix(); },
        "generating active_dof_idx_mat index matrix");

    // Assemble direct history and collector-activated conditions from the current model
    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(step_period); },
        "constructing load matrix (node x 6)");

    // Free-body load balancing requires the complete active support state to be empty
    const bool supports_active =
        !model->_data->conditions.get(bc::SUPPORT).empty();

    if (inertia_relief) {
        logging::error(!supports_active,
            "InertiaRelief: requires no active structural supports. "
            "Clear direct SUPPORT history and remove active support collector definitions.");

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
        logging::error(!supports_active,
            "Rebalancing Loads: requires no active structural supports. "
            "Clear direct SUPPORT history and deselect non-empty support collectors.");

        Timer::measure(
            [&]() { fem::rebalance_loads(*model->_data, global_load_mat); },
            "rebalancing of loads");
    }

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat); },
        "building constraints");

    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // Every linear step solves for du, using either the reset reference
    // configuration or the last converged total displacement as its base.
    model::Field reference_displacement = model->_data->create_field_(
        "REFERENCE_DISPLACEMENT", model::FieldDomain::NODE, 6, false);
    reference_displacement.set_zero();

    // Always linearize around the prepared base configuration. RESUME=NO
    // restored reference positions before run(); RESUME=YES preserved u0.
    const auto& current = *model->_data->positions;
    const auto& reference = *model->_data->positions_reference;
    bool has_deformed_base = false;
    for (Index node = 0; node < reference_displacement.rows; ++node) {
        for (Index dof = 0; dof < 3; ++dof) {
            reference_displacement(node, dof) = current(node, dof) - reference(node, dof);
            has_deformed_base |= reference_displacement(node, dof) != Precision(0);
        }
        for (Index dof = 3; dof < 6; ++dof) {
            reference_displacement(node, dof) = current(node, dof);
            has_deformed_base |= reference_displacement(node, dof) != Precision(0);
        }
    }

    // Constraints act on total displacement: C(u0 + du) = d.
    // Their incremental right-hand side must therefore be d - C u0.
    const Precision tolerance = constraint::NullSpaceOptions{}.feasibility_tolerance;
    for (auto& equation : equations) {
        Precision scale = std::abs(equation.rhs);
        for (const auto& entry : equation.entries) {
            const Precision term =
                entry.coeff * reference_displacement(entry.node_id, entry.dof);
            equation.rhs -= term;
            scale += std::abs(term);
        }
        if (std::abs(equation.rhs) <= tolerance * scale) {
            equation.rhs = Precision(0);
        }
    }

    model::Field reference_internal = model->_data->create_field_(
        "REFERENCE_INTERNAL_FORCES", model::FieldDomain::NODE, 6, false);
    reference_internal.set_zero();

    // K_T is evaluated at the base (u0, T0), while the internal-force anchor
    // uses the target temperature. The difference of internal forces at fixed
    // u0 supplies only the thermal contribution to the load increment.
    auto K = Timer::measure(
        [&]() {
            return model->build_tangent_stiffness_matrix(
                active_dof_idx_mat,
                reference_internal,
                reference_displacement,
                target_temperature,
                reference_displacement,
                model->_data->temperature.get(),
                nullptr,
                false);
        },
        "constructing base tangent stiffness K_T and internal force");

    auto f_external = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing external load matrix -> active vector"
    );

    auto f_reference = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, reference_internal); },
        "reducing reference internal force -> active vector"
    );

    // The tangent assembly evaluated internal force at the target temperature.
    // Compare against the exact base internal force when temperature changes,
    // for both general linear steps and perturbations.
    DynamicVector f_base_internal = f_reference;
    {
        const auto* base_temperature = model->_data->temperature.get();
        bool temperature_changed = (target_temperature == nullptr) != (base_temperature == nullptr);

        if (!temperature_changed && target_temperature) {
            for (Index node = 0; node < target_temperature->rows; ++node) {
                if ((*target_temperature)(node, 0) != (*base_temperature)(node, 0)) {
                    temperature_changed = true;
                    break;
                }
            }
        }

        if (temperature_changed) {
            model::Field base_internal = model->_data->create_field_(
                "BASE_INTERNAL_FORCES", model::FieldDomain::NODE, 6, false);
            base_internal.set_zero();
            model->build_internal_force(
                base_internal,
                reference_displacement,
                base_temperature,
                &reference_displacement,
                base_temperature,
                false,
                nullptr);
            f_base_internal = mattools::reduce_mat_to_vec(active_dof_idx_mat, base_internal);
        }
    }

    // The actual last committed loading is independent of mutable definitions
    // and collector activation. Independent steps instead start at zero.
    DynamicVector f_start = DynamicVector::Zero(f_external.size());
    if (resume && step_start_load) {
        f_start = mattools::reduce_mat_to_vec(active_dof_idx_mat, *step_start_load);
    }

    // Strict increment: no inherited equilibrium residual is corrected.
    // Identical loads, temperatures and constraints therefore yield du=0.
    // The difference of internal forces at fixed u0 is the thermal increment.
    DynamicVector f = (f_external - f_start) - (f_reference - f_base_internal);

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

    // A deformed/prestressed base may have an indefinite or nonsymmetric tangent.
    const auto direct_matrix_type =
        perturbation || has_deformed_base || constraint_method == ConstraintTransformer::Method::Lagrange
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

    // Reject nonfinite system coefficients before passing them to the solver.
    logging::error(mattools::all_finite(A),
        "Matrix A contains NaN/Inf entries");
    logging::error(b.allFinite(),
        "b contains NaN/Inf entries");

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b, direct_matrix_type); },
        "solving constraint system");

    auto u = Timer::measure(
        [&]() { return transformer->recover_displacement(q); },
        "recovering full displacement vector u");

    // The tangent predicts only the internal-force increment. The total
    // internal force must be continued from the last committed field rather
    // than evaluated anew at the deformed base, which would insert the
    // nonlinear residual of earlier linear steps.
    DynamicVector starting_internal = f_base_internal;
    if (resume && start_static_response) {
        const auto& previous = start_static_response->fields[
            static_cast<std::size_t>(io::writer::OutputField::INTERNAL_FORCES)];
        if (previous) starting_internal =
            mattools::reduce_mat_to_vec(active_dof_idx_mat, *previous);
    }

    auto internal_active = Timer::measure(
        [&]() {
            const DynamicVector increment_internal = f_reference - f_base_internal + K * u;
            return perturbation ? DynamicVector(increment_internal)
                                : DynamicVector(starting_internal + increment_internal);
        },
        "computing total or incremental internal nodal forces");

    auto r_support = Timer::measure(
        [&]() { return transformer->support_reactions(K, f, q); },
        "computing support reactions via multipliers (C_supp^T lambda)");

    if (perturbation) {
        // The increment system already supplies incremental reactions and
        // internal forces. Only the external-force output needs conversion.
        if (step_start_load) global_load_mat -= *step_start_load;
    } else {
        // Reactions are continued consistently with the committed internal
        // force, not with a fresh nonlinear evaluation at the base geometry.
        r_support += starting_internal - f_start;
    }

    auto increment_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement increment to matrix form");
    model::Field total_disp_mat = increment_disp_mat;
    total_disp_mat += reference_displacement;

    auto global_react_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r_support); },
        "expanding support reactions to matrix form");

    auto global_internal_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, internal_active); },
        "expanding internal nodal forces to matrix form");

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

    auto reaction_masked = mattools::mask_field(
        global_react_mat,
        support_mask,
        "REACTION_FORCES"
    );

    Timer::measure(
        [&]() {
            using io::writer::OutputField;

            writer->add_loadcase(id, io::writer::WriterStepType::Static);

            // The static step exposes only quantities obtained directly from
            // the solved equilibrium state. Stress, strain, shell resultants,
            // section forces and shear flow are recovered lazily by the output
            // handler if the active request set actually needs them.
            if (perturbation) {
                output.set_perturbation_base(&reference_displacement);
            } else {
                output.set_linearization_base(&reference_displacement);
                output.set_static_response_base(start_static_response.get());
            }
            output.begin_frame();
            output.provide(OutputField::DISPLACEMENT,
                perturbation ? increment_disp_mat : total_disp_mat);
            output.provide(OutputField::EXTERNAL_FORCES,     global_load_mat);
            output.provide(OutputField::INTERNAL_FORCES,     global_internal_mat);
            output.provide(OutputField::REACTION_FORCES,     reaction_masked);
            if (target_temperature) {
                output.provide(OutputField::TEMPERATURE,
                    *target_temperature);
            }

            output.write_frame(*writer, model->_data.get());
        },
        "resolving and writing requested result fields");

    transformer->post_check_static(K, f, q);

    // The constitutive totals must be saved even if the output request did
    // not ask for them. Future resumed steps may request them for the first
    // time, and cannot reconstruct the intervening affine history from u alone.
    if (!perturbation) {
        auto response = std::make_shared<io::writer::StaticResponse>(
            output.capture_static_response());
        response->fields[static_cast<std::size_t>(io::writer::OutputField::INTERNAL_FORCES)] =
            global_internal_mat;
        response->fields[static_cast<std::size_t>(io::writer::OutputField::REACTION_FORCES)] =
            reaction_masked;

        commit_state(total_disp_mat);
        committed_structural_load = std::make_shared<model::Field>(global_load_mat);
        committed_static_response = std::move(response);
    }
    model->step_end();
}

} // namespace loadcase
} // namespace fem
