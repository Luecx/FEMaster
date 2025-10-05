/**
 * @file support.cpp
 * @brief Implements the support boundary condition logic.
 *
 * Support objects traverse their target regions and emit constraint equations
 * that later enforce prescribed displacements or rotations in the solver.
 *
 * @see src/bc/support.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "support.h"

#include "../data/node_data_dict.h"
#include "../model/element/element.h"
#include "../model/model_data.h"

#include <cmath>

namespace fem {
namespace bc {

/**
 * @copydoc Support::Support(NodeRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : node_region(std::move(node_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

/**
 * @copydoc Support::Support(ElementRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : element_region(std::move(element_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

/**
 * @copydoc Support::Support(SurfaceRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : surface_region(std::move(surface_region)), values(values), coordinate_system(std::move(coordinate_system)) {}

/**
 * @copydoc Support::apply
 */
void Support::apply(model::ModelData& model_data, constraint::Equations& equations) {
    if (node_region) {
        for (ID node_id : *node_region) {
            apply_to_node(model_data, equations, node_id);
        }
    } else if (element_region) {
        for (ID element_id : *element_region) {
            if (!model_data.elements[element_id]) {
                continue;
            }
            auto& element = model_data.elements[element_id];
            for (ID node_id : *element) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    } else if (surface_region) {
        for (ID surface_id : *surface_region) {
            if (!model_data.surfaces[surface_id]) {
                continue;
            }
            auto surface = model_data.surfaces[surface_id];
            for (ID node_id : *surface) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    }
}

/**
 * @copydoc Support::apply_to_node
 */
void Support::apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id) {
    Vec6 position_vec = model_data.get(model::POSITION).row(node_id);
    Vec3 position = position_vec.head<3>();

    for (int i = 0; i < 6; ++i) {
        if (std::isnan(values[i])) {
            continue;
        }

        if (coordinate_system) {
            auto local_position = coordinate_system->to_local(position);
            auto transformation = coordinate_system->get_axes(local_position);
            Vec3 direction = transformation.col(i % 3);

            constraint::EquationEntry entry_x = {node_id, static_cast<Dim>((i / 3) * 3), direction[0]};
            constraint::EquationEntry entry_y = {node_id, static_cast<Dim>((i / 3) * 3 + 1), direction[1]};
            constraint::EquationEntry entry_z = {node_id, static_cast<Dim>((i / 3) * 3 + 2), direction[2]};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry_x, entry_y, entry_z});
        } else {
            constraint::EquationEntry entry = {node_id, static_cast<Dim>(i), 1.0};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry}, values[i]);
        }
    }
}

} // namespace bc
} // namespace fem
