/**
 * @file model_thermal_load.cpp
 * @brief Implements thermal boundary-condition and Robin-equation assembly.
 */

#include "model.h"

#include "../bc/neumann/neumann.h"
#include "../bc/robin/robin.h"
#include "../core/logging.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::model {

ThermalLoadData Model::build_thermal_loads(
    const std::vector<std::string>& load_sets,
    Precision                       time
) {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    Field rhs{
        "THERMAL_LOAD",
        FieldDomain::NODE,
        _data->field_rows(FieldDomain::NODE),
        1
    };
    rhs.set_zero();

    bc::RobinEquations equations{};

    for (const std::string& name : load_sets) {
        logging::error(_data->load_cols.has(name),
            "Model: thermal load collector ", name, " does not exist");

        const auto collector = _data->load_cols.get(name);
        logging::error(collector != nullptr,
            "Model: thermal load collector ", name, " is not initialized");

        for (const auto& load : collector->entries()) {
            if (!load) continue;

            if (auto neumann = std::dynamic_pointer_cast<bc::ThermalNeumann>(load)) {
                neumann->apply(*_data, rhs, time, false);
                continue;
            }

            if (auto robin = std::dynamic_pointer_cast<bc::Robin>(load)) {
                robin->apply(*_data, rhs, equations, time, false);
                continue;
            }

            logging::error(false,
                "Model: thermal load collector ", name,
                " contains non-thermal condition ", load->str());
        }
    }

    return {std::move(rhs), std::move(equations)};
}

/**
 * Builds the complete linear thermal operator from element conductivity and
 * additive symbolic Robin equation rows.
 */
SparseMatrix Model::build_thermal_matrix(
    SystemDofIds&              indices,
    const bc::RobinEquations&  equations
) {
    logging::error(indices.cols() == 1,
        "Model: thermal system mapping must have exactly one component");

    SparseMatrix matrix = build_conductivity_matrix(indices);
    TripletList robin_terms{};

    for (const auto& equation : equations) {
        logging::error(equation.row_node_id >= 0
                    && static_cast<Eigen::Index>(equation.row_node_id) < indices.rows(),
            "Model: Robin equation row references an invalid node");
        logging::error(static_cast<Eigen::Index>(equation.row_dof) < indices.cols(),
            "Model: Robin equation row references an invalid DOF");

        const int row = indices(equation.row_node_id, equation.row_dof);
        logging::error(row >= 0,
            "Model: Robin equation targets a thermally inactive node ",
            equation.row_node_id);

        for (const auto& entry : equation.entries) {
            logging::error(entry.node_id >= 0
                        && static_cast<Eigen::Index>(entry.node_id) < indices.rows(),
                "Model: Robin equation references an invalid node");
            logging::error(static_cast<Eigen::Index>(entry.dof) < indices.cols(),
                "Model: Robin equation references an invalid DOF");

            const int col = indices(entry.node_id, entry.dof);
            logging::error(col >= 0,
                "Model: Robin equation references thermally inactive node ",
                entry.node_id);

            robin_terms.emplace_back(row, col, entry.coeff);
        }
    }

    if (!robin_terms.empty()) {
        SparseMatrix robin(matrix.rows(), matrix.cols());
        robin.setFromTriplets(robin_terms.begin(), robin_terms.end());
        matrix += robin;
        matrix.makeCompressed();
    }

    return matrix;
}

} // namespace fem::model
