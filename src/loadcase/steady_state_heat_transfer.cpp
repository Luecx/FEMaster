/**
 * @file steady_state_heat_transfer.cpp
 * @brief Solves linear stationary heat transfer on the compiled FEM model.
 *
 * The formulation integrates heat conduction in the reference configuration
 * and adds Newton-cooling surface contributions to the same scalar operator.
 * Absolute nodal temperatures are prescribed through affine constraint
 * equations. The existing sparse solver and constraint transformer solve the
 * resulting linear system; no separate thermal linear algebra is introduced.
 *
 * @see SteadyStateHeatTransfer
 * @see model::Model::build_thermal_load_matrix
 */

#include "steady_state_heat_transfer.h"

#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/all_finite.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"

#include <limits>
#include <memory>
#include <utility>

namespace fem::loadcase {

using constraint::ConstraintTransformer;

/**
 * Solves the end-of-step stationary balance (K_cond + K_film) T = q.
 *
 * One temperature DOF is numbered for each node participating in a thermal
 * element. Surface heat flux supplies prescribed input, and convection adds
 * both a Robin matrix and an ambient-temperature source. Constraints prescribe
 * absolute nodal temperatures, so the recovered solution is a total field.
 * The result is placed in ModelData::temperature_target for the parser to commit
 * on successful completion, and for later structural steps to reuse.
 *
 * Only the end of the step is evaluated: Abaqus steady-state increment controls
 * are accepted by the reader but do not cause intermediate thermal solves.
 */
void SteadyStateHeatTransfer::run() {
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "STEADY-STATE HEAT TRANSFER");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    logging::error(!perturbation,
        "STEADYHEATTRANSFER does not support perturbation steps");
    logging::error(constraint_method != Method::Lagrange || method != solver::INDIRECT,
        "LAGRANGE constraints require a DIRECT solver");

    // Bind sections before integrating conductivity and surface operators
    model->assign_sections();
    model->step_begin();

    // One active scalar equation per node connected to a thermal element
    const SystemDofIds dofs = Timer::measure(
        [&]() { return model->build_thermal_dof_index_matrix(); },
        "numbering thermal degrees of freedom");

    const Index n_dofs = dofs.size() == 0 ? 0 : static_cast<Index>(dofs.maxCoeff() + 1);
    logging::error(n_dofs > 0,
        "STEADYHEATTRANSFER requires at least one thermal element");

    // Assemble the symmetric reference conduction and convection operators
    auto K = Timer::measure(
        [&]() { return model->build_thermal_conductivity_matrix(dofs); },
        "assembling thermal conductivity");
    K += Timer::measure(
        [&]() { return model->build_thermal_boundary_matrix(dofs, target_conditions, step_period); },
        "assembling convection operator");
    K.makeCompressed();

    // Heat input and ambient convection sources enter the thermal right-hand side
    auto loads = Timer::measure(
        [&]() { return model->build_thermal_load_matrix(target_conditions, step_period, Precision(1)); },
        "assembling thermal loads");
    const DynamicVector f = mattools::reduce_mat_to_vec(dofs, loads);

    // Prescribed absolute temperatures constrain the scalar thermal unknowns
    const auto equations = Timer::measure(
        [&]() { return model->collect_thermal_constraints(target_conditions); },
        "collecting prescribed temperatures");

    ConstraintTransformer::Options options;
    options.method = constraint_method;
    ConstraintTransformer transformer(equations, dofs, n_dofs, options);

    SparseMatrix       A = transformer.assemble_system_matrix(K);
    const DynamicVector b = transformer.assemble_system_rhs(K, f);

    logging::error(mattools::all_finite(A),
        "Thermal system matrix contains NaN/Inf entries");
    logging::error(b.allFinite(),
        "Thermal right-hand side contains NaN/Inf entries");

    // Solve the constrained equilibrium and recover absolute nodal temperatures
    const auto matrix_type = constraint_method == Method::Lagrange
        ? solver::DirectSolverMatrixType::General
        : solver::DirectSolverMatrixType::SPD;
    const DynamicVector q = Timer::measure(
        [&]() { return solve(device, method, A, b, matrix_type); },
        "solving steady-state heat transfer");
    const DynamicVector T = transformer.recover_displacement(q);

    model::Field temperature = mattools::expand_vec_to_mat(dofs, T);
    temperature.name = "TEMPERATURE";

    // Preserve unavailable temperatures on nodes with no thermal equation
    for (Index node = 0; node < temperature.rows; ++node) {
        if (dofs(node, 0) < 0) {
            temperature(node, 0) = std::numeric_limits<Precision>::quiet_NaN();
        }
    }

    // Publish the total temperature for output and subsequent structural steps
    model->_data->temperature_target = std::make_shared<model::Field>(std::move(temperature));

    // HFL is recovered on demand from the newly solved absolute temperature
    using io::writer::OutputField;
    writer->add_loadcase(id, io::writer::WriterStepType::Static);
    output.begin_frame();
    output.provide(OutputField::TEMPERATURE, *model->_data->temperature_target);
    output.write_frame(*writer, model->_data.get());

    transformer.post_check_static(K, f, q);
    model->step_end();
}

} // namespace fem::loadcase
