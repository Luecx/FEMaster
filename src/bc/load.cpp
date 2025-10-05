/**
 * @file load.cpp
 * @brief Implements the concrete load boundary conditions.
 *
 * Each load type converts high-level user input into equivalent nodal loads
 * or boundary-condition entries by traversing the associated FEM regions.
 *
 * @see src/bc/load.h
 * @see src/bc/load_collector.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load.h"

#include "../model/element/element.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <cmath>

namespace fem {
namespace bc {

/**
 * @copydoc CLoad::apply
 */
void CLoad::apply(model::ModelData& model_data, NodeData& bc) {
    (void)model_data;
    for (auto& node_id : *region) {
        for (Dim i = 0; i < 6; ++i) {
            if (!std::isnan(values[i])) {
                bc(node_id, i) += values[i];
            }
        }
    }
}

/**
 * @copydoc DLoad::apply
 */
void DLoad::apply(model::ModelData& model_data, NodeData& bc) {
    auto& node_positions = model_data.get(model::POSITION);
    for (auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_dload(node_positions, bc, values);
    }
}

/**
 * @copydoc PLoad::apply
 */
void PLoad::apply(model::ModelData& model_data, NodeData& bc) {
    auto& node_positions = model_data.get(model::POSITION);
    for (auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_pload(node_positions, bc, pressure);
    }
}

/**
 * @copydoc VLoad::apply
 */
void VLoad::apply(model::ModelData& model_data, NodeData& bc) {
    for (auto& el_id : *region) {
        auto& el_ptr = model_data.elements[el_id];
        auto structural = el_ptr->as<model::StructuralElement>();
        if (structural) {
            structural->apply_vload(bc, values);
        }
    }
}

/**
 * @copydoc TLoad::apply
 */
void TLoad::apply(model::ModelData& model_data, NodeData& bc) {
    for (auto& element_ptr : model_data.elements) {
        if (auto structural = element_ptr->as<model::StructuralElement>()) {
            structural->apply_tload(bc, *temp_field, ref_temp);
        }
    }
}

} // namespace bc
} // namespace fem
