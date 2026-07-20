/**
 * @file support.cpp
 * @brief Implements region expansion and equation generation for supports.
 *
 * A support resolves its active node, element or surface target to individual
 * node identifiers and emits one constraint equation for every prescribed
 * generalized degree of freedom. Global supports map directly to single global
 * DOFs, while oriented supports express a selected local axis as a linear
 * combination of the three corresponding global components.
 *
 * @see support.h
 * @see ../constraints/types/equation.h
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
 * Constructs a support targeting an explicit node region.
 *
 * Prescribed generalized values are applied directly to every node. An optional
 * coordinate system resolves local axes at each node during equation generation.
 *
 * @param node_region Nodes receiving the constraints.
 * @param values Prescribed values for six generalized degrees of freedom.
 * @param coordinate_system Optional coordinate system for oriented constraints.
 */
Support::Support(NodeRegionPtr node_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : node_region_(std::move(node_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Constructs a support targeting an element region.
 *
 * Connected nodes are resolved when the support is applied, preserving the
 * element-based input semantics without duplicating topology data.
 *
 * @param element_region Elements whose connected nodes are constrained.
 * @param values Prescribed values for six generalized degrees of freedom.
 * @param coordinate_system Optional coordinate system for oriented constraints.
 */
Support::Support(ElementRegionPtr element_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : element_region_(std::move(element_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Constructs a support targeting a surface region.
 *
 * Connected nodes are resolved during equation generation, after which the
 * surface-based definition is emitted as node-level constraints.
 *
 * @param surface_region Surfaces whose connected nodes are constrained.
 * @param values Prescribed values for six generalized degrees of freedom.
 * @param coordinate_system Optional coordinate system for oriented constraints.
 */
Support::Support(SurfaceRegionPtr surface_region, const Vec6& values, cos::CoordinateSystem::Ptr coordinate_system)
    : surface_region_(std::move(surface_region)),
      values_(values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Expands the support target and generates node-level constraint equations.
 *
 * Node regions are iterated directly; element and surface regions are reduced
 * to their connected nodes before every node is passed to apply_to_node.
 *
 * @param model_data Model topology and nodal position data.
 * @param equations Equation collection receiving the constraints.
 */
void Support::apply(model::ModelData& model_data, constraint::Equations& equations) {
    if (node_region_) {
        // A node region already contains the exact entities on which the
        // kinematic prescription acts, so no topology lookup is necessary.
        for (ID node_id : *node_region_) {
            apply_to_node(model_data, equations, node_id);
        }
    } else if (element_region_) {
        // An element-region support acts on every node connected to every
        // selected element. Nodes shared by multiple elements are encountered
        // repeatedly because this function deliberately performs no deduplication.
        for (ID element_id : *element_region_) {
            // Ignore empty element slots instead of dereferencing a null model
            // entry during incomplete model construction.
            if (!model_data.elements[element_id]) {
                continue;
            }

            // Iterate the element's connectivity through its common interface
            // and generate the active equations for each connected node.
            auto& element = model_data.elements[element_id];
            for (ID node_id : *element) {
                apply_to_node(model_data, equations, node_id);
            }
        }
    } else if (surface_region_) {
        // A surface-region support expands through surface connectivity in the
        // same manner as the element path.
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
 * Generates the active support equations for one node.
 *
 * Each finite prescribed component becomes a constraint equation. For oriented
 * supports, the coordinate-system basis at the node maps local components into
 * the corresponding global degrees of freedom.
 *
 * @param model_data Nodal position field used for oriented axes.
 * @param equations Equation collection receiving the constraints.
 * @param node_id Global node identifier to constrain.
 */
void Support::apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id) {
    // Position is required to evaluate the axes of a potentially
    // position-dependent coordinate system.
    logging::error(model_data.positions != nullptr, "positions field not set in model data");

    // The position field exposes six-component rows for generalized nodal data;
    // only the first three translational entries form the geometric point.
    const Vec6 position_vec = model_data.positions->row_vec6(static_cast<Index>(node_id));
    const Vec3 position     = position_vec.head<3>();

    // Traverse translations first and rotations second using the common
    // generalized ordering `[Ux, Uy, Uz, Rx, Ry, Rz]`.
    for (int i = 0; i < 6; ++i) {
        // `NaN` is the sentinel for a free degree of freedom and therefore does
        // not generate any equation.
        if (std::isnan(values_[i])) {
            continue;
        }

        if (coordinate_system_) {
            // Convert the global node position into the coordinate system's
            // local parameterization before evaluating its basis matrix.
            const Vec3 local_position = coordinate_system_->to_local(position);
            const auto transformation = coordinate_system_->get_axes(local_position);

            // `i % 3` selects the requested local x, y or z axis. The same axes
            // are used for both generalized blocks; `i / 3` selects translation
            // or rotation in the global equation entries.
            const Vec3 direction      = transformation.col(i % 3);

            // Express one local component as the dot product of the selected
            // axis with the three corresponding global components. The equation
            // constructor used here receives no explicit right-hand side, so the
            // oriented branch creates the default homogeneous constraint; the
            // finite value in `values_[i]` currently acts as the activation flag.
            constraint::EquationEntry entry_x = {node_id, static_cast<Dim>((i / 3) * 3), direction[0]};
            constraint::EquationEntry entry_y = {node_id, static_cast<Dim>((i / 3) * 3 + 1), direction[1]};
            constraint::EquationEntry entry_z = {node_id, static_cast<Dim>((i / 3) * 3 + 2), direction[2]};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry_x, entry_y, entry_z}, values_[i]);
        } else {
            // A global support constrains exactly one generalized DOF. Its
            // coefficient is one and the prescribed value becomes the equation
            // right-hand side.
            constraint::EquationEntry entry = {node_id, static_cast<Dim>(i), 1.0};
            equations.emplace_back(std::initializer_list<constraint::EquationEntry>{entry}, values_[i]);
        }
    }
}

/**
 * Builds the diagnostic representation of the support definition.
 *
 * The result identifies the active node, element or surface region, reports its
 * cardinality and lists the six nominal prescribed values.
 *
 * @return Human-readable support description.
 */
std::string Support::str() const {
    // Resolve the one active region pointer into a type-qualified name and
    // cardinality. The lambda keeps this selection local to string generation.
    const auto region_name = [&]() -> std::string {
        if (node_region_)
            return "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        if (surface_region_)
            return "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        if (element_region_)
            return "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        return "(unknown)";
    }();

    // Labels follow the storage order in `values_` and make the printed
    // prescription independent of numeric DOF indices.
    static constexpr const char* labels[] = {
        "Ux", "Uy", "Uz", "Rx", "Ry", "Rz"
    };

    std::ostringstream os;
    os << "Support: target=" << region_name << ", dof=[";

    // Append only finite prescribed entries and manage separators without a
    // trailing comma.
    bool first = true;
    for (int i = 0; i < 6; ++i) {
        if (std::isnan(values_[i]))
            continue;

        if (!first)
            os << ", ";

        os << labels[i] << '=' << values_[i];
        first = false;
    }

    os << ']';

    // Identify the local basis only when the prescription is not expressed
    // directly in global coordinates.
    if (coordinate_system_)
        os << ", orientation=" << coordinate_system_->name;

    return os.str();
}
} // namespace bc
} // namespace fem
