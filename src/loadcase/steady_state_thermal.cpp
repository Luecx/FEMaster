/**
 * @file steady_state_thermal.cpp
 * @brief Implements linear steady-state heat-conduction analysis.
 *
 * The analysis assembles the stationary thermal balance
 *
 *     (K_T + K_h) T = f_q + f_h,
 *
 * where K_T is element conductivity, K_h is the convection boundary matrix,
 * f_q is prescribed surface heat input and f_h is the ambient convection source.
 * Prescribed absolute temperatures are enforced by the common constraint
 * transformer. The converged temperature field is subsequently differentiated
 * by each thermal element to recover conductive heat flux.
 *
 * @see SteadyStateThermal
 * @see model::Model::build_conductivity_matrix
 * @see model::Model::compute_heat_flux
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "steady_state_thermal.h"

#include "../bc/neumann/convection.h"
#include "../bc/neumann/heat_flux.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../model/model.h"

#include <cmath>
#include <memory>

namespace fem::loadcase {

using constraint::ConstraintTransformer;

/**
 * Solves the constrained stationary thermal system and writes its result fields.
 *
 * The processing sequence is:
 *
 * 1. Assign sections and enumerate one scalar temperature DOF per active node.
 * 2. Collect prescribed temperatures from the selected support collectors.
 * 3. Assemble element conductivity and thermal surface conditions.
 * 4. Transform and solve the affine constrained linear system.
 * 5. Recover nodal temperature and integration-point conductive heat flux.
 * 6. Write one static thermal result frame.
 *
 * Thermal load collectors must contain only `HeatFlux` and `Convection`
 * conditions. This explicit selection prevents mechanical loads from being
 * interpreted as scalar heat input. Prescribed temperatures may be absent when
 * convection already anchors the absolute temperature level; otherwise the
 * solver reports the singular unanchored conductivity system.
 */
void SteadyStateThermal::run() {
    // Print the analysis header and validate parser-assigned dependencies.
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

    // Bind material-bearing sections before evaluating thermal elements.
    model->assign_sections();

    // Enumerate one scalar unknown at every thermally active node.
    auto thermal_dof_ids = Timer::measure(
        [&]() { return model->build_thermal_index_matrix(); },
        "generating thermal DOF index matrix"
    );
    logging::error(thermal_dof_ids.rows() > 0,
        "STEADYSTATETHERMAL: compiled model contains no nodes");
    logging::error(thermal_dof_ids.maxCoeff() >= 0,
        "STEADYSTATETHERMAL: model contains no thermally active element nodes");

    const Index        active_temperatures = static_cast<Index>(thermal_dof_ids.maxCoeff() + 1);
    const Eigen::Index system_size         = static_cast<Eigen::Index>(active_temperatures);

    // Collect only prescribed-temperature entries from common support collectors.
    auto groups = Timer::measure(
        [&]() { return model->collect_temperature_constraints(supps); },
        "collecting prescribed temperatures"
    );
    report_constraint_groups(groups);

    auto equations = groups.flatten();

    // Assemble element conductivity in the global scalar temperature system.
    auto conductivity = Timer::measure(
        [&]() { return model->build_conductivity_matrix(thermal_dof_ids); },
        "constructing thermal conductivity matrix K_T"
    );
    logging::error(conductivity.rows() == system_size
                && conductivity.cols() == system_size,
        "STEADYSTATETHERMAL: conductivity matrix dimensions do not match thermal DOFs");

    // Assemble scalar nodal heat input and convection boundary-matrix triplets.
    model::Field nodal_heat_source{
        "THERMAL_LOAD",
        model::FieldDomain::NODE,
        model->_data->field_rows(model::FieldDomain::NODE),
        1
    };
    nodal_heat_source.set_zero();

    TripletList boundary_terms{};
    Timer::measure(
        [&]() {
            for (const std::string& name : loads) {
                logging::error(model->_data->load_cols.has(name),
                    "STEADYSTATETHERMAL: load collector ", name, " does not exist");

                const auto collector = model->_data->load_cols.get(name);
                logging::error(collector != nullptr,
                    "STEADYSTATETHERMAL: load collector ", name, " is not initialized");

                for (const auto& load : collector->entries()) {
                    if (!load) continue;

                    const bool thermal = std::dynamic_pointer_cast<bc::HeatFlux>(load) != nullptr
                                      || std::dynamic_pointer_cast<bc::Convection>(load) != nullptr;
                    logging::error(thermal,
                        "STEADYSTATETHERMAL: load collector ", name,
                        " contains non-thermal condition ", load->str());

                    load->apply(
                        *model->_data,
                        nodal_heat_source,
                        Precision(0),
                        false,
                        &thermal_dof_ids,
                        &boundary_terms
                    );
                }
            }
        },
        "assembling thermal boundary conditions"
    );

    // Convert nodal scalar storage to the compact active thermal vector.
    DynamicVector heat_source = DynamicVector::Zero(system_size);
    for (Index node = 0; node < nodal_heat_source.rows; ++node) {
        const int thermal_dof = thermal_dof_ids(node, 0);
        if (thermal_dof >= 0) {
            heat_source(thermal_dof) += nodal_heat_source(node, 0);
        }
    }

    // Add the convection matrix to element conductivity after duplicate surface
    // contributions have been accumulated by sparse triplet construction.
    if (!boundary_terms.empty()) {
        SparseMatrix boundary_matrix(system_size, system_size);
        boundary_matrix.setFromTriplets(boundary_terms.begin(), boundary_terms.end());
        conductivity += boundary_matrix;
        conductivity.makeCompressed();
    }

    // Validate solver and constraint-backend compatibility.
    if (constraint_method == ConstraintTransformer::Method::Lagrange && method == solver::INDIRECT) {
        logging::error(false,
            "STEADYSTATETHERMAL: LAGRANGE constraints require the DIRECT solver");
    }

    const auto direct_matrix_type =
        constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

    // Construct the affine representation of prescribed absolute temperatures.
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

    // Assemble the transformed matrix and right-hand side, including the
    // inhomogeneous prescribed-temperature contribution.
    auto system_matrix = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(conductivity); },
        "assembling constrained thermal matrix"
    );
    auto system_rhs = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(conductivity, heat_source); },
        "assembling constrained thermal RHS"
    );

    // Reject invalid element or transformation results before solver entry.
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

    // A fully prescribed reduced system has no remaining algebraic unknowns and
    // therefore requires no sparse solver call.
    auto solution = Timer::measure(
        [&]() -> DynamicVector {
            if (system_matrix.rows() == 0) {
                return DynamicVector{};
            }
            return solve(device, method, system_matrix, system_rhs, direct_matrix_type);
        },
        "solving steady-state thermal system"
    );

    // Recover the complete active temperature vector from transformed coordinates.
    const DynamicVector active_temperature = transformer->recover_displacement(solution);
    logging::error(active_temperature.size() == system_size,
        "STEADYSTATETHERMAL: recovered temperature vector has the wrong size");
    logging::error(active_temperature.allFinite(),
        "STEADYSTATETHERMAL: recovered temperature contains NaN or Inf");

    // Expand active temperatures into an explicit scalar nodal result field.
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

    // Apply Fourier's law to the converged temperature field in each thermal element.
    auto heat_flux = Timer::measure(
        [&]() { return model->compute_heat_flux(temperature); },
        "computing conductive heat flux"
    );

    // Write temperature to native and FRD output. Integration-point heat flux is
    // retained by the native result format and ignored by the nodal-only FRD writer.
    Timer::measure(
        [&]() {
            writer->add_loadcase(id, io::writer::WriterStepType::Static);
            writer->write_field(temperature, "TEMPERATURE", model->_data.get());
            writer->write_field(heat_flux  , "HEAT_FLUX" , model->_data.get());
        },
        "writing thermal result fields"
    );

    transformer->post_check_static(conductivity, heat_source, solution);
}

} // namespace fem::loadcase
