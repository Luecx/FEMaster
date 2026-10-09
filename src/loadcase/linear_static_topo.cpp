/**
 * @file linear_static_topo.cpp
 * @brief Implements topology-aware linear static analysis.
 *
 * Topology density and orientation settings affect structural stiffness and
 * sensitivity recovery. Current ModelData conditions supply loads and prescribed
 * components using the same assembly path as LinearStatic. Temporary scaling
 * fields belong to this numerical solve; direct history and named definitions
 * remain model-owned. Constraint reduction and sparse solution reuse existing
 * structural analysis utilities.
 *
 * @see LinearStaticTopo
 * @see LinearStatic
 */

#include "linear_static_topo.h"

#include "../constraints/transformer/constraint_transformer.h"
#include "../constraints/types/equation.h"
#include "../constraints/types/rbm.h"
#include "../core/logging.h"
#include "../mattools/assemble.h"
#include "../mattools/mask_field.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../solve/eigval/solve_eigval.h"
#include "tools/inertia_relief.h"
#include "tools/rebalance_loads.h"

#include <cmath>
#include <limits>

using fem::constraint::ConstraintTransformer;

namespace fem { namespace loadcase {

/**
 * Executes density/orientation-dependent equilibrium for current model conditions.
 *
 * Topology parameters establish element stiffness scaling and material bases.
 * Current direct conditions and collector-activated conditions supply the structural RHS
 * and prescribed components. Free-body balancing validates actual support state
 * and may temporarily add an RBM for inertia relief during constraint collection.
 *
 * The reduced linear solve uses the scaled stiffness operator. Global recovery
 * supplies displacement, stress, strain and topology sensitivities for requested
 * output. Temporary model field handles and element caches are restored after
 * writing; persistent condition history remains owned by ModelData.
 */
void LinearStaticTopo::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "LINEAR STATIC TOPO");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Prepare section assignments and analysis-local element state
    model->assign_sections();

    logging::error(density != nullptr,
        "LinearStaticTopo: density field not initialized");

    model::Field::Ptr stiffness = std::make_shared<model::Field>(
        "TOPO_DENSITY_STIFFNESS", model::FieldDomain::ELEMENT, density->rows, 1);

    for (Index i = 0; i < stiffness->rows; ++i) {
        if (!std::isnan((*density)(i)) && std::isfinite((*density)(i))) {
            (*stiffness)(i) = std::pow((*density)(i), exponent);
        }
    }

    model->_data->element_stiffness_scale = stiffness;
    model->_data->material_orientation    = orientation;
    model->step_begin();

    // Enumerate the current structural system in global generalized coordinates
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return model->build_structural_dof_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // Assemble direct history and collector-activated conditions from the current model
    auto global_load_mat = Timer::measure(
        [&]() { return model->build_load_matrix(); },
        "constructing load matrix (node x 6)"
    );

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

                // This RBM is analysis-local bookkeeping, not a Model factory
                // concern. It is stored directly and removed after collection.
                model->_data->rbms.emplace_back(model->_data->elem_sets.get("EALL"));
            },
            "InertiaRelief: adjusting external load matrix and adding RBM"
        );
    }

    if (rebalance_loads) {
        logging::error(!supports_active,
            "Rebalancing Loads: requires no active structural supports. "
            "Clear direct SUPPORT history and deselect non-empty support collectors.");

        Timer::measure(
            [&]() { fem::rebalance_loads(*model->_data, global_load_mat); },
            "rebalancing of loads"
        );
    }

    auto groups = Timer::measure(
        [&]() { return model->collect_constraints(active_dof_idx_mat); },
        "building constraints"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    auto K = Timer::measure(
        [&]() { return model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing stiffness matrix K(rho^p, theta)"
    );

    model::Field reference_displacement{
        "REFERENCE_DISPLACEMENT",
        model::FieldDomain::NODE,
        model->_data->field_rows(model::FieldDomain::NODE),
        6
    };
    reference_displacement.set_zero();

    model::Field reference_internal{
        "REFERENCE_INTERNAL_FORCES",
        model::FieldDomain::NODE,
        model->_data->field_rows(model::FieldDomain::NODE),
        6
    };
    model->build_internal_force(
        reference_internal,
        reference_displacement,
        model->_data->temperature.get(),
        nullptr,
        nullptr,
        false,
        nullptr);

    auto f_external = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load matrix -> active external vector"
    );
    auto f_reference = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, reference_internal); },
        "reducing reference internal force -> active vector"
    );
    auto f = f_external - f_reference;

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

    auto CT = Timer::measure(
        [&]() {
            ConstraintTransformer::Options copt;
            copt.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations,
                active_dof_idx_mat,
                K.rows(),
                copt
            );
        },
        "building constraint transformer"
    );

    logging::info(true, "");
    logging::info(true, "Constraint summary");
    logging::up();
    logging::info(true, "m (rows of C)        : ", CT->report().equations);
    logging::info(true, "n (cols of C)        : ", CT->report().dofs);
    if (CT->rank_known()) {
        logging::info(true, "rank(C)              : ", CT->rank());
    } else {
        logging::info(true, "rank(C)              : not computed");
    }
    logging::info(true, "method               : ", CT->method_name());
    logging::info(true, "solver unknowns      : ", CT->unknowns());
    logging::info(true, "homogeneous          : ", CT->homogeneous() ? "true" : "false");
    logging::info(true, "feasible             : ", CT->feasible() ? "true" : "false");
    logging::info(!CT->feasible(), "residual ||C u - d|| : ", CT->report().residual_norm);
    logging::down();

    if (inertia_relief) {
        logging::error(!model->_data->rbms.empty(),
            "InertiaRelief: expected a temporary RBM to be present, but rbms is empty");
        model->_data->rbms.pop_back();
    }

    auto A = Timer::measure(
        [&]() { return CT->assemble_system_matrix(K); },
        "assembling constraint system matrix"
    );
    auto b = Timer::measure(
        [&]() { return CT->assemble_system_rhs(K, f); },
        "assembling constraint system RHS"
    );

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
        logging::error(!badA,
            "Matrix A contains NaN/Inf entries");
        logging::error(b.allFinite(),
            "b contains NaN/Inf entries");
    }

    auto q = Timer::measure(
        [&]() { return solve(device, method, A, b, direct_matrix_type); },
        "solving reduced system A q = b"
    );

    auto u = Timer::measure(
        [&]() { return CT->recover_displacement(q); },
        "recovering full displacement vector u"
    );
    auto r_support = Timer::measure(
        [&]() { return CT->support_reactions(K, f, q); },
        "computing support reactions via multipliers (C_supp^T lambda)"
    );

    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, u); },
        "expanding displacement vector to matrix form"
    );
    auto global_react_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, r_support); },
        "expanding support reactions to matrix form"
    );

    // Compliance is also required to form the topology density sensitivity,
    // so compute it once here and expose it as an already available field.
    // Other ordinary FE results remain lazy output-handler recoveries.
    auto compliance_raw = model->compute_compliance(global_disp_mat);
    model::Field density_grad = model->_data->create_field_(
        "DENS_GRAD", model::FieldDomain::ELEMENT, 1, false);

    const Index element_count = static_cast<Index>(model->_data->elements.size());
    for (Index r = 0; r < element_count; ++r) {
        const Precision rho  = (*density)(r, 0);
        density_grad(r, 0) = -exponent * compliance_raw(r, 0) / rho;
    }

    const bool has_orientation = orientation != nullptr;

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

    using io::writer::OutputField;

    // The historical topology defaults include orientation fields only when an
    // orientation field is configured. Re-applying defaults here is safe:
    // explicit input-deck requests are retained by OutputRequestHandler.
    if (has_orientation) {
        output.set_defaults({
            OutputField::DISPLACEMENT,
            OutputField::STRAIN,
            OutputField::STRESS,
            OutputField::EXTERNAL_FORCES,
            OutputField::REACTION_FORCES,
            OutputField::COMPLIANCE,
            OutputField::DENS_GRAD,
            OutputField::VOLUME,
            OutputField::DENSITY,
            OutputField::SHEAR_FLOW,
            OutputField::ORIENTATION_GRAD,
            OutputField::ORIENTATION
        });
    }

    writer->add_loadcase(id, io::writer::WriterStepType::Static);

    output.begin_frame();
    output.provide(OutputField::DISPLACEMENT,        global_disp_mat);
    output.provide(OutputField::EXTERNAL_FORCES,     global_load_mat);
    output.provide(OutputField::REACTION_FORCES,     reaction_masked);
    output.provide(OutputField::COMPLIANCE,          compliance_raw);
    output.provide(OutputField::DENS_GRAD,           density_grad);
    output.provide(OutputField::DENSITY,             *density);
    if (has_orientation) output.provide(OutputField::ORIENTATION, *orientation);

    output.write_frame(*writer, model->_data.get());

    CT->post_check_static(K, f, q);

    model->_data->element_stiffness_scale = nullptr;
    model->_data->material_orientation    = nullptr;
    model->step_end();
}

}} // namespace fem::loadcase
