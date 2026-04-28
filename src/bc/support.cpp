/**
 * @file support.cpp
 * @brief Implements the support boundary condition logic.
 *
 * Support objects traverse their target regions and emit constraint equations
 * that later enforce prescribed displacements or rotations in the solver.
 *
 * @see src/bc/support.h
 * @see src/constraints/types/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "support.h"

#include "../data/field.h"
#include "../model/element/element.h"
#include "../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem {
namespace bc {
/**
 * @copydoc Support::Support(NodeRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : node_region_(std::move(node_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * @copydoc Support::Support(ElementRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : element_region_(std::move(element_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * @copydoc Support::Support(SurfaceRegionPtr,const Vec6&,cos::CoordinateSystem::Ptr)
 */
Support::Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : surface_region_(std::move(surface_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * @copydoc Support::apply
 */
void Support::apply(model::ModelData& model_data, constraint::Equations& equations) {
    if (node_region_) {
        // Node supports can be emitted directly because the target region
        // already contains the final constrained node ids.
        for (ID node_id : *node_region_) {
            apply_to_node(model_data, equations, node_id);
        }
    } else if (element_region_) {
        // Element supports constrain every node of every selected element.
        // Missing element pointers are skipped so partially initialized models
        // do not crash during reporting or intermediate setup.
        for (ID element_id : *element_region_) {
            if (!model_data.elements[element_id]) {
                continue;
            }

            auto& element = model_data.elements[element_id];
            for (ID node_id : *element) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    } else if (surface_region_) {
        // Surface supports work like element supports, but traverse surface
        // connectivity instead of element connectivity.
        for (ID surface_id : *surface_region_) {
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
    logging::error(model_data.positions != nullptr, "positions field not set in model data");

    const Vec6 position_vec = model_data.positions->row_vec6(static_cast<Index>(node_id));
    const Vec3 position     = position_vec.head<3>();

    for (int i = 0; i < 6; ++i) {
        if (std::isnan(values_[i])) {
            continue;
        }

        if (coordinate_system_) {
            // Local supports constrain the global components that represent
            // the requested local axis at this node position.
            const Vec3 local_position = coordinate_system_->to_local(position);
            const auto transformation = coordinate_system_->get_axes(local_position);
            const Vec3 direction      = transformation.col(i % 3);

            constraint::EquationEntry entry_x = {node_id, static_cast<Dim>((i / 3) * 3), direction[0]};
            constraint::EquationEntry entry_y = {node_id, static_cast<Dim>((i / 3) * 3 + 1), direction[1]};
            constraint::EquationEntry entry_z = {node_id, static_cast<Dim>((i / 3) * 3 + 2), direction[2]};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry_x, entry_y, entry_z});
        } else {
            // Global supports map one prescribed DOF directly to one equation.
            constraint::EquationEntry entry = {node_id, static_cast<Dim>(i), 1.0};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry}, values_[i]);
        }
    }
}

std::string Support::str() const {
    auto region_name = [&]() -> std::string {
        if (node_region_) {
            return std::string("NSET ") + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        }

        if (surface_region_) {
            return std::string("SFSET ") + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        }

        if (element_region_) {
            return std::string("ELSET ") + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        }

        return std::string("(unknown)");
    }();

    static const char* labels[6] = {"Ux", "Uy", "Uz", "Rx", "Ry", "Rz"};

    std::ostringstream os;
    os << "Support: target="
       << region_name
       << ", dof=[";

    bool first = true;
    for (int i = 0; i < 6; ++i) {
        if (!std::isnan(values_[i])) {
            if (!first) {
                os << ", ";
            }

            os << labels[i]
               << "="
               << values_[i];
            first = false;
        }
    }

    if (first) {
        os << "(none)";
    }

    if (coordinate_system_) {
        os << ", orientation=" << coordinate_system_->name;
    }

    return os.str();
}
} // namespace bc
} // namespace fem
