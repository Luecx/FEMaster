/**
 * @file steady_state_thermal.cpp
 * @brief Implements linear steady-state heat-conduction analysis.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "steady_state_thermal.h"

#include "../bc/neumann/thermal_load.h"
#include "../core/logging.h"
#include "../core/timer.h"
#include "../model/element/element_thermal.h"
#include "../model/model.h"

#include <cmath>
#include <memory>

namespace fem::loadcase {

using constraint::ConstraintTransformer;

namespace {

/**
 * Projects element integration-point heat flux to a writer-compatible nodal field.
 *
 * The current thermal element interface recovers HFL at integration points. For
 * visualization and FRD output, one mean flux vector is first assigned to every
 * node of the contributing element. The established element-nodal averaging
 * routine then combines connected-element contributions at each global node.
 */
model::Field project_heat_flux_to_nodes(model::ModelData& data,
                                        const model::Field& integration_flux) {
    logging::error(integration_flux.domain == model::FieldDomain::ELEMENT_IP,
        "HFL projection requires an ELEMENT_IP source field");
    logging::error(integration_flux.components >= 3,
        "HFL projection requires three heat-flux components");

    const Index element_count = static_cast<Index>(data.elements.size());
    model::Field element_flux{
        "HFL_ELEMENT_NODAL",
        model::FieldDomain::ELEMENT_NODAL,
        data.field_rows(model::FieldDomain::ELEMENT_NODAL),
        3
    };
    model::Field element_weights{
        "HFL_ELEMENT_WEIGHTS",
        model::FieldDomain::ELEMENT,
        element_count,
        1
    };
    element_flux.set_zero();
    element_weights.set_zero();

    for (const auto& element : data.elements) {
        if (!element || !element->as<model::ThermalElement>()) continue;

        const Index ip_count = element->num_ip();
        if (ip_count == 0) continue;

        Vec3 mean_flux = Vec3::Zero();
        for (Index ip = 0; ip < ip_count; ++ip) {
            const Index row = element->ip_index(ip);
            for (Dim component = 0; component < 3; ++component) {
                mean_flux(component) += integration_flux(row, component);
            }
        }
        mean_flux /= static_cast<Precision>(ip_count);

        for (Index local_node = 0; local_node < element->n_nodes(); ++local_node) {
            const Index row = static_cast<Index>(element->elem_nodal_offset) + local_node;
            for (Dim component = 0; component < 3; ++component) {
                element_flux(row, component) = mean_flux(component);
            }
        }
        element_weights(static_cast<Index>(element->elem_id), 0) = Precision(1);
    }

    auto nodal_flux = data.element_nodal_to_nodal(element_flux, element_weights, "HFL");
    nodal_flux.check_finite("Nodal heat flux");
    return nodal_flux;
}

} // namespace

void SteadyStateThermal::run() {
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

    model->assign_sections();

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

    auto groups = Timer::measure(
        [&]() { return model->collect_temperature_constraints(supps); },
        "collecting prescribed temperatures"
    );
    report_constraint_groups(groups);
    auto equations = groups.flatten();
    logging::error(!equations.empty(),
        "STEADYSTATETHERMAL: at least one prescribed temperature is required");

    // Conductive element contribution.
    auto conductivity = Timer::measure(
        [&]() { return model->build_conductivity_matrix(thermal_dof_ids); },
        "constructing thermal conductivity matrix K_T"
    );
    logging::error(conductivity.rows() == active_temperatures
                && conductivity.cols() == active_temperatures,
        "STEADYSTATETHERMAL: conductivity matrix dimensions do not match thermal DOFs");

    // Natural thermal boundary conditions. Prescribed heat flux contributes only
    // to f_T; convection contributes both K_h and f_h.
    DynamicVector heat_source = DynamicVector::Zero(active_temperatures);
    TripletList boundary_terms;

    Timer::measure(
        [&]() {
            for (const auto& name : loads) {
                logging::error(model->_data->thermal_load_cols.has(name),
                    "STEADYSTATETHERMAL: thermal load collector ", name, " does not exist");
                const auto collector = model->_data->thermal_load_cols.get(name);
                logging::error(collector != nullptr,
                    "STEADYSTATETHERMAL: thermal load collector ", name, " is not initialized");

                for (const auto& load : collector->entries()) {
                    if (load) {
                        load->apply(*model->_data, thermal_dof_ids, boundary_terms, heat_source);
                    }
                }
            }
        },
        "assembling thermal Neumann conditions"
    );

    if (!boundary_terms.empty()) {
        SparseMatrix boundary_matrix(active_temperatures, active_temperatures);
        boundary_matrix.setFromTriplets(boundary_terms.begin(), boundary_terms.end());
        conductivity += boundary_matrix;
        conductivity.makeCompressed();
    }

    if (constraint_method == ConstraintTransformer::Method::Lagrange && method == solver::INDIRECT) {
        logging::error(false,
            "STEADYSTATETHERMAL: LAGRANGE constraints require the DIRECT solver");
    }

    const auto direct_matrix_type =
        constraint_method == ConstraintTransformer::Method::Lagrange
            ? solver::DirectSolverMatrixType::General
            : solver::DirectSolverMatrixType::SPD;

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

    auto system_matrix = Timer::measure(
        [&]() { return transformer->assemble_system_matrix(conductivity); },
        "assembling constrained thermal matrix"
    );
    auto system_rhs = Timer::measure(
        [&]() { return transformer->assemble_system_rhs(conductivity, heat_source); },
        "assembling constrained thermal RHS"
    );

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

    auto solution = Timer::measure(
        [&]() {
            if (system_matrix.rows() == 0) {
                return DynamicVector(0);
            }
            return solve(device, method, system_matrix, system_rhs, direct_matrix_type);
        },
        "solving steady-state thermal system"
    );

    const DynamicVector active_temperature = transformer->recover_displacement(solution);
    logging::error(active_temperature.size() == active_temperatures,
        "STEADYSTATETHERMAL: recovered temperature vector has the wrong size");
    logging::error(active_temperature.allFinite(),
        "STEADYSTATETHERMAL: recovered temperature contains NaN or Inf");

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

    // Element HFL is currently recovered at integration points. Convert it to a
    // NODE field before handing it to writers, which also addresses FRD output.
    auto integration_flux = Timer::measure(
        [&]() { return model->compute_heat_flux(temperature); },
        "computing conductive heat flux"
    );
    auto heat_flux = Timer::measure(
        [&]() { return project_heat_flux_to_nodes(*model->_data, integration_flux); },
        "projecting heat flux to nodes"
    );

    // Reaction heat flow is the full thermal residual. Free nodes are zero up to
    // solver tolerance; prescribed-temperature nodes retain the required power.
    const DynamicVector active_reaction = conductivity * active_temperature - heat_source;
    logging::error(active_reaction.allFinite(),
        "STEADYSTATETHERMAL: recovered reaction heat flux contains NaN or Inf");

    model::Field reaction{
        "RFL",
        model::FieldDomain::NODE,
        model->_data->field_rows(model::FieldDomain::NODE),
        1
    };
    reaction.set_zero();
    for (Index node = 0; node < reaction.rows; ++node) {
        const int thermal_dof = thermal_dof_ids(node, 0);
        if (thermal_dof >= 0) {
            reaction(node, 0) = active_reaction(thermal_dof);
        }
    }

    Timer::measure(
        [&]() {
            writer->add_loadcase(id, io::writer::WriterStepType::Static);
            writer->write_field(temperature, "TEMPERATURE", model->_data.get());
            writer->write_field(heat_flux  , "HFL"        , model->_data.get());
            writer->write_field(reaction   , "RFL"        , model->_data.get());
        },
        "writing thermal result fields"
    );

    transformer->post_check_static(conductivity, heat_source, solution);
}

} // namespace fem::loadcase
