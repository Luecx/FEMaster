/**
 * @file linear_static.cpp
 * @brief Implements the linear static load case leveraging constraint maps.
 *
 * The analysis assembles step-owned loads and supports, then reduces
 * K_T du = delta f through the chosen constraint map and recovers
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
 * Ordinary linear statics solves total equilibrium about the reference
 * state. Perturbations alone linearize about a committed, prestressed base.
 *
 * The constraint map solves total linear equilibrium for general steps and
 * tangent increments for perturbations with the same configured backend.
 * Expanded global displacements drive stress/strain and reaction recovery, then
 * requested result fields are written and analysis-local element caches released.
 * The analysis evaluates its step-owned initial and target conditions without
 * mutating the physical definitions.
 */
void LinearStatic::run() {
    // Apply step-specific defaults only when no output variables were requested
    if (output.empty()) {
        using io::writer::OutputField;
        output.request({
            {OutputField::DISPLACEMENT        , model::FieldDomain::NODE},
            {OutputField::STRAIN              , model::FieldDomain::NODE},
            {OutputField::STRESS              , model::FieldDomain::NODE},
            {OutputField::STRESS_TOP          , model::FieldDomain::NODE},
            {OutputField::STRESS_BOT          , model::FieldDomain::NODE},
            {OutputField::SHELL_RESULTANTS    , model::FieldDomain::NODE},
            {OutputField::EXTERNAL_FORCES     , model::FieldDomain::NODE},
            {OutputField::REACTION_FORCES     , model::FieldDomain::NODE},
            {OutputField::LOCAL_SECTION_FORCES, model::FieldDomain::ELEMENT_NODAL},
            {OutputField::SHEAR_FLOW          , model::FieldDomain::ELEMENT_NODAL}
        });
    }

    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Ordinary linear statics never continues a deformed tangent: the
    // reference configuration is its geometry throughout assembly.
    // Perturbations retain the prestressed configuration unmodified.
    if (!perturbation) {
        *model->_data->positions = *model->_data->positions_reference;
    }

    // error checks if we are in a perturbation analysis
    logging::error(!perturbation || model->_data->contacts.empty(),
        "LINEARSTATIC PERTURBATION: contact tangent is not supported");
    logging::error(!perturbation || (!inertia_relief && !rebalance_loads),
        "LINEARSTATIC PERTURBATION: inertia relief and rebalancing are not supported");
    logging::error(!perturbation || method == solver::DIRECT,
        "LINEARSTATIC PERTURBATION: a DIRECT solver is required for general tangents");
    logging::error(subcases.empty() || perturbation,
        "LINEARSTATIC: SUBCASE requires a PERTURBATION step");
    logging::error(constraint_method != Method::Lagrange || method != solver::INDIRECT,
        "LAGRANGE constraints require a DIRECT solver");
    logging::error(!inertia_relief || (model->_data->elem_sets.has("EALL") && model->_data->elem_sets.get("EALL") != nullptr),
        "InertiaRelief: EALL element set is not available");
    logging::error(!inertia_relief || target_conditions.get(bc::SUPPORT).empty(),
       "InertiaRelief: requires no active structural supports. "
       "Clear direct SUPPORT history and remove active support collector definitions.");
    logging::error(!rebalance_loads || target_conditions.get(bc::SUPPORT).empty(),
        "Rebalancing Loads: requires no active structural supports. "
        "Clear direct SUPPORT history and deselect non-empty support collectors.");

    model->assign_sections();
    model->step_begin();

    // Structural temperatures are prescribed fields, not thermal constraints.
    // The solver evaluates thermal forces at the step's target temperature.
    model::Field* target_temperature = model->_data->temperature_target
        ? model->_data->temperature_target.get()
        : model->_data->temperature.get();

    // Enumerate the current structural system in global generalized coordinates
    auto system_dofs = Timer::measure(
        [&]() { return model->build_structural_dofs(); },
        "numbering structural degrees of freedom");

    // Assemble step initial/target conditions and selected collectors. Each
    // subcase evaluates its own load target against the common initial state.
    auto load_mats = Timer::measure(
        [&]() {
            std::vector<model::Field> mats;
            if (subcases.empty()) {
                mats.push_back(build_step_load_matrix(step_period, Precision(1), model::SystemDofs{}, nullptr));
            }
            for (const auto& [name, conditions] : subcases) {
                mats.push_back(model->build_load_matrix(
                    initial_conditions, conditions, step_period, Precision(1), model::SystemDofs{}, nullptr));
            }
            return mats;
        },
        "constructing load matrix (node x 6)");

    // Inertia relief and rebalancing are rejected for perturbations and
    // therefore only ever act on the single non-subcase load field.
    auto& global_load_mat = load_mats.front();

    if (inertia_relief) {
        // Inertia relief adds one temporary RBM directly to the same
        // collection used by parser-created constraints. It is removed
        // again immediately after constraint collection below.
        Timer::measure(
            [&]() {fem::apply_inertia_relief(*model->_data, global_load_mat, inertia_relief_consider_point_masses);
                model->_data->rbms.emplace_back(model->_data->elem_sets.get("EALL"));
            },"InertiaRelief: adjusting external load matrix and adding RBM");
    }

    if (rebalance_loads) {
        Timer::measure(
            [&]() { fem::rebalance_loads(*model->_data, global_load_mat); },
            "rebalancing of loads");
    }

    // collector all the constraints that are active in this analysis
    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(system_dofs, target_conditions); },
        "building constraints");

    // report if required by the user and turn into equations
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // The total displacement is solved from the reference in ordinary statics.
    // Only perturbations use the committed configuration as tangent base.
    model::Field base_displacement = model->_data->create_field_(
        "LINEARIZATION_DISPLACEMENT", model::FieldDomain::NODE, 6, false);
    base_displacement.set_zero();

    if (perturbation) {
        const auto& current   = *model->_data->positions;
        const auto& reference = *model->_data->positions_reference;
        for (Index node = 0; node < base_displacement.rows; ++node) {
            for (Index dof = 0; dof < 3; ++dof) {
                base_displacement(node, dof) = current(node, dof) - reference(node, dof);
            }
            for (Index dof = 3; dof < 6; ++dof) {
                base_displacement(node, dof) = current(node, dof);
            }
        }

        // Incremental boundary condition: C du = d - C u0.
        const Precision tolerance = constraint::NullSpaceOptions{}.feasibility_tolerance;
        for (auto& equation : equations) {
            Precision scale = std::abs(equation.rhs);
            for (const auto& entry : equation.entries) {
                const Precision term =
                    entry.coeff * base_displacement(entry.node_id, entry.dof);
                equation.rhs -= term;
                scale += std::abs(term);
            }
            if (std::abs(equation.rhs) <= tolerance * scale) equation.rhs = Precision(0);
        }
    }

    model::Field reference_internal = model->_data->create_field_(
        "REFERENCE_INTERNAL_FORCES", model::FieldDomain::ELEMENT_NODAL, 6, false);

    // Tangent at (0,Tref) for ordinary statics or (u0,T0) for perturbation.
    // The internal-force anchor uses the target temperature.
    const model::Field* base_temperature = perturbation
        ? model->_data->temperature.get()
        : model->_data->temperature_reference.get();

    auto K = Timer::measure(
        [&]() {
            return model->build_tangent_stiffness_matrix(
                system_dofs, reference_internal, base_displacement,
                target_temperature, base_displacement, base_temperature,
                nullptr, false);
        },
        "constructing tangent stiffness and reference internal force");

    model::Field reference_nodal = reference_internal.nodal_sum(*model->_data);
    model->build_contact_force(system_dofs, reference_nodal, &K);

    const Index n_rhs = static_cast<Index>(load_mats.size());

    auto f_reference = Timer::measure(
        [&]() { return system_dofs.reduce(reference_nodal); },
        "reducing reference internal force -> active vector");

    // Perturbations linearize about the committed base. Its internal force
    // and the committed start load are shared by every right-hand side.
    DynamicVector f_base_internal = f_reference;
    DynamicVector f_start         = DynamicVector::Zero(f_reference.size());
    if (perturbation) {
        if (step_start_load) {
            f_start = system_dofs.reduce(*step_start_load);
        }

        model::Field base_internal = model->_data->create_field_(
            "BASE_INTERNAL_FORCES", model::FieldDomain::ELEMENT_NODAL, 6, false);
        model->build_internal_force(
            base_internal, base_displacement, base_temperature,
            &base_displacement, base_temperature, false);
        f_base_internal =
            system_dofs.reduce(base_internal.nodal_sum(*model->_data));
    }

    // General linear statics: K0 u = F - fint(0,Ttarget).
    // Perturbation: K(u0,T0) du = dF - [fint(u0,Ttarget)-fint(u0,T0)].
    DynamicMatrix f = Timer::measure(
        [&]() {
            DynamicMatrix rhs(f_reference.size(), n_rhs);
            for (Index i = 0; i < n_rhs; ++i) {
                const DynamicVector f_external = system_dofs.reduce(load_mats[i]);
                rhs.col(i) = perturbation
                    ? DynamicVector((f_external - f_start) - (f_reference - f_base_internal))
                    : DynamicVector(f_external - f_reference);
            }
            return rhs;
        },
        "reducing external load matrices -> active vectors");

    // Convert external forces to increments for perturbation output
    if (perturbation && step_start_load) {
        for (auto& load_mat : load_mats) load_mat -= *step_start_load;
    }

    // Only perturbations may inherit a prestressed/indefinite tangent.
    const auto direct_matrix_type =
        perturbation || constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations,
                system_dofs,
                K.rows(),
                options);
        },
        "building constraint transformer");

    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)        : ", transformer->report().equations);
    logging::info(true, "n (cols of C)        : ", transformer->report().dofs);
    logging::info(true, "rank(C)              : ", transformer->rank_known () ? std::to_string(transformer->rank()) : "not computed");
    logging::info(true, "method               : ", transformer->method_name());
    logging::info(true, "solver unknowns      : ", transformer->unknowns   ());
    logging::info(true, "homogeneous          : ", transformer->homogeneous() ? "true" : "false");
    logging::info(true, "feasible             : ", transformer->feasible   () ? "true" : "false");
    logging::info(true, "residual ||C u - d|| : ", transformer->feasible   () ? "within tolerance" : std::to_string(transformer->report().residual_norm));
    logging::info(true, "right-hand sides     : ", n_rhs);
    logging::down();

    // Remove the temporary inertia-relief RBM after constraint collection
    if (inertia_relief)
        model->_data->rbms.pop_back();

    auto A = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(K); },
        "assembling constraint system matrix");

    auto b = Timer::measure(
        [&]() {
            DynamicMatrix rhs(transformer->unknowns(), n_rhs);
            for (Index i = 0; i < n_rhs; ++i) {
                rhs.col(i) = transformer->assemble_system_rhs(K, DynamicVector(f.col(i)));
            }
            return rhs;
        },
        "assembling constraint system RHS");

    // Reject nonfinite system coefficients before passing them to the solver.
    logging::error(mattools::all_finite(A),
        "Matrix A contains NaN/Inf entries");
    logging::error(b.allFinite(),
        "b contains NaN/Inf entries");

    // One factorization serves every subcase right-hand side.
    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b, direct_matrix_type); },
        "solving constraint system");

    // export stiffness files in .mtx format
    if (!stiffness_file.empty()) {
        io::writer::write_mtx      (stiffness_file + "_K.mtx", K);
        io::writer::write_mtx      (stiffness_file + "_A.mtx", A);
        io::writer::write_mtx_dense(stiffness_file + "_b.mtx", b);
    }

    const auto support_mask = groups.support_mask(system_dofs);

    writer->add_loadcase(id, io::writer::WriterStepType::Static);

    // Recover results consistently with the linearization used for K.
    if (perturbation) {
        output.set_perturbation_base(&base_displacement);
    } else {
        output.set_linearization_base(&base_displacement, base_temperature);
    }

    for (Index i = 0; i < n_rhs; ++i) {
        const DynamicVector f_i = f.col(i);
        const DynamicVector q_i = q.col(i);

        if (!subcases.empty()) {
            logging::info(true, "");
            logging::info(true, "Subcase ", i + 1, ": ", subcases[i].first);
        }

        auto u = Timer::measure(
            [&]() { return transformer->recover_displacement(q_i); },
            "recovering full displacement vector u");

        // The solver's affine internal forces are the only force outputs.
        // Perturbations emit their increment, not the prestressed base force.
        auto internal_active = Timer::measure(
            [&]() {
                return perturbation
                    ? DynamicVector(f_reference - f_base_internal + K * u)
                    : DynamicVector(f_reference + K * u);
            },
            "computing affine internal nodal forces");

        auto r_support = Timer::measure(
            [&]() { return transformer->support_reactions(K, f_i, q_i); },
            "computing support reactions");

        auto global_disp_mat = Timer::measure(
            [&]() { return system_dofs.expand(u); },
            "expanding solved displacement to nodal field");

        auto global_react_mat = Timer::measure(
            [&]() { return system_dofs.expand(r_support); },
            "expanding support reactions to matrix form");

        auto global_internal_mat = Timer::measure(
            [&]() { return system_dofs.expand(internal_active); },
            "expanding internal nodal forces to matrix form");

        auto reaction_masked = mattools::mask_field(
            global_react_mat, support_mask, "REACTION_FORCES"
        );

        Timer::measure(
            [&]() {
                using io::writer::OutputField;

                // Subcases are written as numbered frames of one load case
                if (subcases.empty()) {
                    output.begin_frame();
                } else {
                    output.begin_frame(Precision(i + 1), "_" + std::to_string(i + 1));
                }
                output.provide(OutputField::DISPLACEMENT   , global_disp_mat);
                output.provide(OutputField::EXTERNAL_FORCES, load_mats[i]);
                output.provide(OutputField::INTERNAL_FORCES, global_internal_mat);
                output.provide(OutputField::REACTION_FORCES, reaction_masked);
                if (target_temperature) {
                    output.provide(OutputField::TEMPERATURE,*target_temperature);
                }

                output.write_frame(*writer, model->_data.get());
            },
            "resolving and writing requested result fields");

        transformer->post_check_static(K, f_i, q_i);

        // Commit the total displacement only for ordinary linear statics,
        // which never carry subcases and therefore solve exactly once.
        if (!perturbation) {
            commit_state(global_disp_mat);
            committed_structural_load = std::make_shared<model::Field>(load_mats[i]);
        }
    }
    model->step_end();
}

} // namespace loadcase
} // namespace fem
