/**
 * @file steady_state_thermal.cpp
 * @brief Implements linear steady-state heat-conduction analysis.
 *
 * The implementation enumerates one temperature unknown per thermally active
 * node, assembles the global conductivity matrix, imposes prescribed absolute
 * temperatures through the common constraint transformer and solves the zero-
 * source stationary balance. The recovered temperature drives element-level
 * Fourier heat-flux evaluation before both result fields are written.
 *
 * Thermal loads and capacity terms are deliberately absent. They belong to
 * future source/flux/convection conditions and transient thermal analysis.
 *
 * @see SteadyStateThermal
 * @see model::Model::build_conductivity_matrix
 * @see model::Model::compute_heat_flux
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#include "steady_state_thermal.h"

#include "../core/logging.h"
#include "../core/timer.h"
#include "../model/model.h"

#include <cmath>
#include <memory>

namespace fem::loadcase {

using constraint::ConstraintTransformer;

/**
 * Solves the constrained stationary heat equation and writes thermal results.
 *
 * The analysis performs the following phases:
 *
 * - assign compiled solid sections and enumerate scalar thermal unknowns,
 * - collect only `bc::Temperature` entries from the selected support collectors,
 * - assemble the conductivity operator and its zero heat-source right-hand side,
 * - transform and solve the affine constrained system,
 * - recover a one-component nodal temperature field,
 * - evaluate integration-point heat flux through Fourier's law,
 * - write `TEMPERATURE` and `HEAT_FLUX` as one static result frame.
 *
 * A missing temperature constraint leaves the zero-source conductivity problem
 * singular and is rejected explicitly. Additional disconnected components still
 * require one prescribed temperature each; the selected sparse solver diagnoses
 * any remaining singularity.
 */
void SteadyStateThermal::run() {
    // Print the analysis header and validate parser-assigned dependencies
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "STEADY-STATE THERMAL ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    logging::error(model != nullptr,
        "STEADYSTATETHERMAL: model is not initialized");
    logging::error(writer != nullptr,
        "STEADYSTATETHERMAL: result writer is not initialized");

    // Bind material-bearing sections before thermal element evaluation
    model->assign_sections();

    // Enumerate one active temperature unknown per thermal element node
    auto thermal_dof_ids = Timer::measure(
        [&]() { return model->build_thermal_index_matrix(); },
        "generating thermal DOF index matrix"
    );
    logging::error(thermal_dof_ids.rows() > 0,
        "STEADYSTATETHERMAL: compiled model contains no nodes");
    const int max_thermal_dof = thermal_dof_ids.maxCoeff();
    logging::error(max_thermal_dof >= 0,
        "STEADYSTATETHERMAL: model contains no thermally active element nodes");
    const Index active_temperatures = static_cast<Index>(max_thermal_dof + 1);

    // Collect absolute temperature constraints from common support collectors
    auto groups = Timer::measure(
        [&]() { return model->collect_temperature_constraints(supps); },
        "collecting prescribed temperatures"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();
    logging::error(!equations.empty(),
        "STEADYSTATETHERMAL: at least one prescribed temperature is required");

    // Assemble the symmetric conduction operator and the zero source vector
    auto conductivity = Timer::measure(
        [&]() { return model->build_conductivity_matrix(thermal_dof_ids); },
        "constructing thermal conductivity matrix K_T"
    );
    logging::error(conductivity.rows() == active_temperatures
                && conductivity.cols() == active_temperatures,
        "STEADYSTATETHERMAL: conductivity matrix dimensions do not match thermal DOFs");

    const DynamicVector heat_source = DynamicVector::Zero(active_temperatures);

    // Validate solver and constraint-backend compatibility
    if (constraint_method == ConstraintTransformer::Method::Lagrange && method == solver::INDIRECT) {
        logging::error(false,
            "STEADYSTATETHERMAL: LAGRANGE constraints require the DIRECT solver");
    }

    const auto direct_matrix_type =
        constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

    // Construct the affine representation of prescribed absolute temperatures
    auto transformer = Timer::measure(
        [&]() {
            ConstraintTransformer::Options options;
            options.method = constraint_method;
            return std::make_unique<ConstraintTransformer>(
                equations,
                thermal_dof_ids,
                active_temperatures,
                options
            );
        },
        "building thermal constraint transformer"
    );
    logging::error(transformer->feasible(),
        "STEADYSTATETHERMAL: prescribed temperatures are inconsistent");

    // Assemble the solver system including affine Dirichlet contributions
    auto system_matrix = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(conductivity); },
        "assembling constrained thermal matrix"
    );
    auto system_rhs = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(conductivity, heat_source); },
        "assembling constrained thermal RHS"
    );

    // Reject invalid element or transformation results before entering a solver
    bool matrix_finite = true;
    for (int column = 0; column < system_matrix.outerSize() && matrix_finite; ++column) {
        for (SparseMatrix::InnerIterator value(system_matrix, column); value; ++value) {
            if (!std::isfinite(value.value())) {
                matrix_finite = false;
                break;
            }
        }
    }
    logging::error(matrix_finite,
        "STEADYSTATETHERMAL: constrained conductivity matrix contains NaN or Inf");
    logging::error(system_rhs.allFinite(),
        "STEADYSTATETHERMAL: constrained right-hand side contains NaN or Inf");

    // Solve the reduced or augmented system. A fully prescribed reduced system
    // has no remaining algebraic unknowns and therefore needs no solver call.
    auto solution = Timer::measure(
        [&]() {
            if (system_matrix.rows() == 0) {
                return DynamicVector(0);
            }
            return solve(device, method, system_matrix, system_rhs, direct_matrix_type);
        },
        "solving steady-state thermal system"
    );

    // Recover the complete active temperature vector from affine coordinates
    const DynamicVector active_temperature = transformer->recover_displacement(solution);
    logging::error(active_temperature.size() == active_temperatures,
        "STEADYSTATETHERMAL: recovered temperature vector has the wrong size");
    logging::error(active_temperature.allFinite(),
        "STEADYSTATETHERMAL: recovered temperature contains NaN or Inf");

    // Expand active values into an explicit scalar nodal result field
    model::Field temperature{
        "TEMPERATURE",
        model::FieldDomain::NODE,
        model->_data->field_rows(model::FieldDomain::NODE),
        1
    };
    temperature.set_zero();
    for (Index node = 0; node < temperature.rows; ++node) {
        const int thermal_dof = thermal_dof_ids(node, 0);
        if (thermal_dof >= 0) {
            temperature(node, 0) = active_temperature(thermal_dof);
        }
    }

    // Recover conductive heat flux from the converged nodal temperatures
    auto heat_flux = Timer::measure(
        [&]() { return model->compute_heat_flux(temperature); },
        "computing conductive heat flux"
    );

    // Write one stationary result frame in native and supported FRD formats
    Timer::measure(
        [&]() {
            writer->add_loadcase(id, io::writer::WriterStepType::Static);
            writer->write_field(temperature, "TEMPERATURE", model->_data.get());
            writer->write_field(heat_flux  , "HEAT_FLUX" , model->_data.get());
        },
        "writing thermal result fields"
    );

    // Verify the recovered solution against the transformed thermal system
    transformer->post_check_static(conductivity, heat_source, solution);
}

} // namespace fem::loadcase
