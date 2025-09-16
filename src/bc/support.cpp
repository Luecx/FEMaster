/******************************************************************************
 * @file support.cpp
 * @brief Implements the Support class for applying boundary conditions in an FEM model.
 *
 * @details The Support class applies constraints to nodes, elements, or surfaces based on
 *          the specified region. The implementation includes handling of transformations
 *          when a local coordinate system is provided.
 *
 * @date Created on 28.11.2024
 * @author Finn Eggers
 ******************************************************************************/

#include "support.h"
#include "../model/model_data.h"
#include "../model/element/element.h"
#include "../data/node_data_dict.h"

namespace fem {

// Constructor for node region
Support::Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : node_region(std::move(node_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

// Constructor for element region
Support::Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : element_region(std::move(element_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

// Constructor for surface region
Support::Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : surface_region(std::move(surface_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

void Support::apply(model::ModelData& model_data, constraint::Equations& equations) {
    if (node_region) {
        for (ID node_id : *node_region) {
            apply_to_node(model_data, equations, node_id);
        }
    } else if (element_region) {
        for (ID el_id : *element_region) {
            if (!model_data.elements[el_id]) continue;
            auto& el = model_data.elements[el_id];
            for(ID node_id: (*el)) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    } else if (surface_region) {
        for (ID el_id : *surface_region) {
            if (!model_data.surfaces[el_id]) continue;
            auto el = model_data.surfaces[el_id];
            for(ID node_id: *el) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    }
}

void Support::apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id) {
    Vec6 position_vec = model_data.get(model::POSITION).row(node_id);
    Vec3 position = position_vec.head(3);

    for (int i = 0; i < 6; i++) {
        if (!std::isnan(values[i])) {
            if (coordinate_system) {
                auto local_position = coordinate_system->to_local(position);
                auto transformation = coordinate_system->get_axes(local_position);
                Vec3 vals = transformation.col(i % 3);

                constraint::EquationEntry en1 = {node_id, (Dim)((i / 3) * 3), vals[0]};
                constraint::EquationEntry en2 = {node_id, (Dim)((i / 3) * 3 + 1), vals[1]};
                constraint::EquationEntry en3 = {node_id, (Dim)((i / 3) * 3 + 2), vals[2]};
                equations.push_back(constraint::Equation({en1, en2, en3}));
            } else {
                constraint::EquationEntry en = {node_id, (Dim)i, 1.0};
                equations.push_back(constraint::Equation({en}, values[i]));
            }
        }
    }
}

} // namespace fem
