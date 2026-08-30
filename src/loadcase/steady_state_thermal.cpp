/**
 * @file steady_state_thermal.cpp
 * @brief Implements linear steady-state heat-conduction analysis.
 *
 * The analysis assembles one scalar temperature degree of freedom per thermally
 * active node. Element conductivity and Robin boundary operators form the global
 * thermal operator, while prescribed heat flux and the ambient part of Robin
 * convection form the thermal right-hand side. Prescribed temperatures are
 * collected as Dirichlet equations and enforced through the common constraint
 * transformer used by the structural solvers.
 *
 * After solving the constrained linear system, the active temperature vector is
 * expanded back to a nodal field. Thermal elements then recover conductive heat
 * flux in global reference coordinates; the element-nodal values are projected
 * to nodes for result output.
 *
 * @see SteadyStateThermal
 * @see model::Model::build_thermal_index_matrix
 * @see model::Model::build_conductivity_matrix
 * @see model::Model::build_thermal_load_matrix
 * @see constraint::ConstraintTransformer
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "steady_state_thermal.h"

#include "../core/logging.h"
#include "../core/timer.h"
#include "../mattools/reduce_mat_to_vec.h"
#include "../model/model.h"

#include <cmath>
#include <memory>

namespace fem::loadcase {

using constraint::ConstraintTransformer;

/**
 * Solves the stationary linear heat-conduction problem for the current loadcase.
 *
 * The solution sequence is deliberately kept explicit and follows the algebraic
 * thermal problem from model data to result fields:
 *
 * - assign compiled sections required by thermal element material access,
 * - enumerate the scalar active temperature unknowns,
 * - collect prescribed-temperature Dirichlet equations,
 * - assemble element conductivity and thermal boundary conditions,
 * - add Robin operator terms to the conductivity matrix,
 * - construct and validate the selected constraint transformation,
 * - solve the constrained linear system,
 * - recover the full nodal temperature field,
 * - recover and project conductive heat flux,
 * - write temperature and heat-flux result fields.
 *
 * The unconstrained steady-state equation is
 *
 * \f[
 *     (K_T + K_R)\,T = f_q + f_R,
 * \f]
 *
 * where `K_T` is the element conductivity matrix, `K_R` contains Robin boundary
 * terms, `f_q` contains prescribed Neumann heat flux and `f_R` contains the
 * prescribed ambient part of Robin conditions. Dirichlet temperatures are not
 * inserted into either side directly; they are represented by affine equations
 * and handled by `ConstraintTransformer`.
 *
 * Conductive heat flux is recovered after convergence from Fourier's law in the
 * reference configuration and is independent of the boundary heat-flow vector.
 */
void SteadyStateThermal::run() {
    // Print the analysis header before any model preparation or timed operation.
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "===============================================================================================");
    logging::info(true, "STEADY-STATE THERMAL ANALYSIS");
    logging::info(true, "===============================================================================================");
    logging::info(true, "");

    // Validate the loadcase resources required by assembly and result output.
    logging::error(model != nullptr,
        "STEADYSTATETHERMAL: model is not initialized");
    logging::error(writer != nullptr,
        "STEADYSTATETHERMAL: result writer is not initialized");

    // Bind compiled elements to their section and material definitions before
    // thermal element matrices query conductivity or other thermal properties.
    model->assign_sections();

    // Enumerate exactly one temperature DOF at every node connected to a thermal
    // element. Inactive nodes remain outside the algebraic thermal system.
    auto thermal_dof_ids = Timer::measure(
        [&]() { return model->build_thermal_index_matrix(); },
        "generating thermal DOF index matrix"
    );
    logging::error(thermal_dof_ids.rows() > 0,
        "STEADYSTATETHERMAL: compiled model contains no nodes");
    logging::error(thermal_dof_ids.cols() == 1,
        "STEADYSTATETHERMAL: thermal DOF index matrix must have one column");
    logging::error(thermal_dof_ids.maxCoeff() >= 0,
        "STEADYSTATETHERMAL: model contains no thermally active element nodes");

    const Index active_temperatures = static_cast<Index>(thermal_dof_ids.maxCoeff() + 1);
    const Eigen::Index system_size  = static_cast<Eigen::Index>(active_temperatures);

    // Collect only prescribed-temperature entries from the selected Dirichlet
    // collectors and flatten them into the common constraint representation.
    auto groups = Timer::measure(
        [&]() { return model->collect_temperature_constraints(supps); },
        "collecting prescribed temperatures"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();

    // Assemble the element conductivity operator over all thermally active
    // elements using the scalar thermal DOF mapping.
    auto conductivity = Timer::measure(
        [&]() { return model->build_conductivity_matrix(thermal_dof_ids); },
        "constructing thermal conductivity matrix K_T"
    );
    logging::error(conductivity.rows() == system_size
                && conductivity.cols() == system_size,
        "STEADYSTATETHERMAL: conductivity matrix dimensions do not match active DOFs");

    // Assemble natural and mixed thermal boundary conditions. Neumann conditions
    // contribute only to the scalar nodal RHS, while Robin conditions contribute
    // both an ambient source term and temperature-dependent operator triplets.
    TripletList boundary_terms{};
    auto nodal_heat_source = Timer::measure(
        [&]() {
            return model->build_thermal_load_matrix(
                loads,
                thermal_dof_ids,
                boundary_terms
            );
        },
        "assembling thermal boundary conditions"
    );

    // Reduce the scalar nodal heat-flow field to the active algebraic ordering.
    auto heat_source = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(thermal_dof_ids, nodal_heat_source); },
        "reducing thermal load field -> active RHS vector"
    );
    logging::error(heat_source.size() == system_size,
        "STEADYSTATETHERMAL: thermal RHS dimensions do not match active DOFs");

    // Add the assembled Robin boundary operator to the element conductivity
    // matrix. Repeated triplets are summed by Eigen during sparse construction.
    if (!boundary_terms.empty()) {
        SparseMatrix boundary_matrix(system_size, system_size);
        boundary_matrix.setFromTriplets(boundary_terms.begin(), boundary_terms.end());
        conductivity += boundary_matrix;
        conductivity.makeCompressed();
    }

    // Lagrange multiplier constraints enlarge the system into an indefinite
    // saddle-point problem and therefore require the direct solver path.
    if (constraint_method == ConstraintTransformer::Method::Lagrange && method == solver::INDIRECT) {
        logging::error(false,
            "STEADYSTATETHERMAL: LAGRANGE constraints require the DIRECT solver");
    }

    // The unconstrained conductivity operator is symmetric positive definite
    // under sufficient thermal restraint. Lagrange augmentation removes the SPD
    // structure and must be announced as a general direct system.
    const auto direct_matrix_type =
        constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

    // Build the affine constraint transformation for the scalar thermal DOFs and
    // reject inconsistent prescribed temperatures before system assembly.
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

    // Transform the full thermal equilibrium equation into the selected
    // constrained representation, including affine Dirichlet contributions.
    auto system_matrix = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(conductivity); },
        "assembling constrained thermal matrix"
    );
    auto system_rhs = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(conductivity, heat_source); },
        "assembling constrained thermal RHS"
    );

    // Validate all assembled coefficients explicitly before invoking a sparse
    // solver so invalid material, geometry or boundary data fail at this level.
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
        "STEADYSTATETHERMAL: constrained thermal matrix contains NaN or Inf");
    logging::error(system_rhs.allFinite(),
        "STEADYSTATETHERMAL: constrained right-hand side contains NaN or Inf");

    // Solve the constrained stationary system. A zero-sized transformed system
    // is valid when all active temperatures are eliminated by constraints.
    auto solution = Timer::measure(
        [&]() -> DynamicVector {
            if (system_matrix.rows() == 0) return DynamicVector{};
            return solve(device, method, system_matrix, system_rhs, direct_matrix_type);
        },
        "solving steady-state thermal system"
    );

    // Recover the complete active temperature vector from the constrained solver
    // coordinates and verify it before expansion into nodal storage.
    const DynamicVector active_temperature = transformer->recover_displacement(solution);
    logging::error(active_temperature.size() == system_size,
        "STEADYSTATETHERMAL: recovered temperature vector has the wrong size");
    logging::error(active_temperature.allFinite(),
        "STEADYSTATETHERMAL: recovered temperature contains NaN or Inf");

    // Expand the algebraic temperature vector into a one-component NODE field.
    // Thermally inactive nodes retain the zero initialization of the expansion.
    auto temperature = Timer::measure(
        [&]() {
            auto field = mattools::expand_vec_to_mat(thermal_dof_ids, active_temperature);
            field.name = "TEMPERATURE";
            return field;
        },
        "expanding thermal solution -> nodal field"
    );

    // Recover Fourier heat flux at each thermal element node and project only
    // thermal-element contributions into a global nodal result field.
    auto heat_flux = Timer::measure(
        [&]() {
            auto element_heat_flux = model->compute_heat_flux(temperature);
            auto mapping_weights   = model->build_field_mapping_weights(false, true);
            return model->_data->element_nodal_to_nodal(
                element_heat_flux,
                mapping_weights,
                "HEAT_FLUX"
            );
        },
        "computing nodal conductive heat flux"
    );

    // Write the converged scalar temperature and global conductive heat-flux
    // fields as one static result step.
    Timer::measure(
        [&]() {
            writer->add_loadcase(id, io::writer::WriterStepType::Static);
            writer->write_field(temperature, "TEMPERATURE", model->_data.get());
            writer->write_field(heat_flux, "HEAT_FLUX", model->_data.get());
        },
        "writing thermal result fields"
    );

    // Evaluate the original full-system equilibrium after solution recovery so
    // constraint enforcement does not hide a thermal residual inconsistency.
    transformer->post_check_static(conductivity, heat_source, solution);
}

} // namespace fem::loadcase
