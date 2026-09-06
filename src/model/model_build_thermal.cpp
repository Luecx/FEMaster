/**
 * @file model_build_thermal.cpp
 * @brief Implements thermal system construction helpers.
 *
 * @author Finn Eggers
 * @date 06.09.2026
 */

#include "../mattools/numerate_dofs.h"
#include "element/element_thermal.h"
#include "model.h"

namespace fem::model {

/**
 * Enumerates the scalar temperature DOFs required by thermal elements.
 *
 * Every node referenced by at least one ThermalElement receives exactly one
 * active temperature component. Nodes referenced only by non-thermal elements
 * remain inactive. The resulting mask is converted to contiguous zero-based
 * system indices.
 *
 * @return Node-by-one matrix of global thermal system DOF identifiers.
 */
SystemDofIds Model::build_thermal_dof_index_matrix() {
    logging::error(_data->positions != nullptr,
        "Model: POSITION field is not initialized");

    SystemDofs mask{_data->positions->rows, 1};
    mask.fill(false);

    for (const auto& element : _data->elements) {
        if (element == nullptr || element->as<ThermalElement>() == nullptr) continue;

        for (ID local_node = 0; local_node < element->n_nodes(); ++local_node) {
            mask(element->nodes()[local_node], 0) = true;
        }
    }

    return mattools::numerate_dofs(mask);
}

} // namespace fem::model
