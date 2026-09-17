/**
 * @file support.cpp
 * @brief Implements target expansion and structural Dirichlet equation generation.
 *
 * A `Support` converts a region-level structural prescription into scalar rows
 * of the global constraint system. Node targets are consumed directly, while
 * element and surface targets are expanded through their compiled connectivity.
 *
 * In the global basis a prescribed component creates a unit row. For a local
 * basis with axis `e_a`, the corresponding translational constraint is
 *
 *     e_a^T u = u_bar,
 *
 * with `u = [Ux, Uy, Uz]^T`; rotational constraints use the same projection on
 * `r = [Rx, Ry, Rz]^T`. Position-dependent coordinate systems are evaluated at
 * each individual node before the row coefficients are formed.
 *
 * @see Support
 * @see Dirichlet
 * @see constraint::Equation
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "support.h"

#include "../../core/logging.h"
#include "../../data/field.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem {
namespace bc {

/**
 * Constructs a support targeting an explicit node region.
 *
 * The supplied generalized values are applied to every node in the region. An
 * optional coordinate system changes only the basis in which those values are
 * interpreted; the final equations are always written in global model DOFs.
 *
 * @param node_region Nodes receiving the structural prescription.
 * @param values Prescribed values ordered as `[Ux, Uy, Uz, Rx, Ry, Rz]` with
 *               `NaN` for free components.
 * @param coordinate_system Optional local basis of the prescribed components.
 */
Support::Support(NodeRegionPtr             node_region,
                 const Vec6&               values,
                 cos::CoordinateSystem::Ptr coordinate_system)
    : node_region_      (std::move(node_region)),
      values_           (values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Constructs a support targeting all nodes connected to an element region.
 *
 * The element definition is retained semantically and expanded only when the
 * support is applied. This avoids copying connectivity into a separate node set
 * and keeps the support consistent with the compiled element topology.
 *
 * @param element_region Elements whose connected nodes receive the prescription.
 * @param values Prescribed values ordered as `[Ux, Uy, Uz, Rx, Ry, Rz]` with
 *               `NaN` for free components.
 * @param coordinate_system Optional local basis of the prescribed components.
 */
Support::Support(ElementRegionPtr          element_region,
                 const Vec6&               values,
                 cos::CoordinateSystem::Ptr coordinate_system)
    : element_region_   (std::move(element_region)),
      values_           (values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Constructs a support targeting all nodes connected to a surface region.
 *
 * The surface connectivity is resolved during application so the stored support
 * remains a surface-level semantic definition while still generating ordinary
 * nodal constraint equations.
 *
 * @param surface_region Surfaces whose connected nodes receive the prescription.
 * @param values Prescribed values ordered as `[Ux, Uy, Uz, Rx, Ry, Rz]` with
 *               `NaN` for free components.
 * @param coordinate_system Optional local basis of the prescribed components.
 */
Support::Support(SurfaceRegionPtr          surface_region,
                 const Vec6&               values,
                 cos::CoordinateSystem::Ptr coordinate_system)
    : surface_region_   (std::move(surface_region)),
      values_           (values),
      coordinate_system_(std::move(coordinate_system)) {}

/**
 * Expands the configured target region and appends its nodal constraint rows.
 *
 * Exactly one semantic target must be active. A node region already contains the
 * final entities. Element and surface regions are expanded through their global
 * connectivity, after which `apply_to_node()` constructs the scalar equations
 * for each prescribed generalized component.
 *
 * This function deliberately performs no node deduplication. If two selected
 * elements or surfaces share a node, both traversals generate the same support
 * rows. Detection or elimination of redundant equations belongs to the global
 * constraint-processing stage rather than to semantic region expansion.
 *
 * @param model_data Compiled model topology and nodal geometry.
 * @param equations Constraint collection receiving the generated rows.
 */
void Support::apply(model::ModelData& model_data, constraint::Equations& equations) {
    // Validate the semantic target before traversing compiled topology
    const int active_regions = static_cast<int>(node_region_    != nullptr)
                             + static_cast<int>(element_region_ != nullptr)
                             + static_cast<int>(surface_region_ != nullptr);

    logging::error(active_regions == 1,
        "SUPPORT: exactly one node, element or surface region must be configured");

    // A direct node target requires no topological expansion
    if (node_region_) {
        for (ID node_id : *node_region_) {
            apply_to_node(model_data, equations, node_id);
        }
        return;
    }

    // Expand every selected element to its connected global nodes. Invalid or
    // uninitialized compiled entries violate the support's target invariant and
    // are reported explicitly rather than being silently ignored.
    if (element_region_) {
        for (ID element_id : *element_region_) {
            logging::error(element_id >= 0 && static_cast<Index>(element_id) < model_data.elements.size(),
                "SUPPORT: element ", element_id, " is outside the compiled element domain");

            auto& element = model_data.elements[static_cast<Index>(element_id)];
            logging::error(element != nullptr,
                "SUPPORT: element ", element_id, " is not initialized");

            for (ID node_id : *element) {
                apply_to_node(model_data, equations, node_id);
            }
        }
        return;
    }

    // Expand every selected surface to its connected global nodes using the same
    // explicit validation of the compiled target
    for (ID surface_id : *surface_region_) {
        logging::error(surface_id >= 0 && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "SUPPORT: surface ", surface_id, " is outside the compiled surface domain");

        auto surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "SUPPORT: surface ", surface_id, " is not initialized");

        for (ID node_id : *surface) {
            apply_to_node(model_data, equations, node_id);
        }
    }
}

/**
 * Converts the six stored prescriptions for one node into scalar constraint rows.
 *
 * For a global component `i`, the generated equation is simply
 *
 *     u_i = value_i.
 *
 * If a coordinate system is assigned, let `A(x)` be the local-to-global basis
 * matrix evaluated at the nodal position and let `e_a = A.col(a)` denote the
 * selected local axis. For a translational component the row becomes
 *
 *     e_a,x Ux + e_a,y Uy + e_a,z Uz = value_i,
 *
 * and for a rotational component the same coefficients multiply
 * `[Rx, Ry, Rz]`. Consequently the equation is expressed entirely in global
 * structural DOFs even though the prescribed scalar component is local.
 *
 * `NaN` values are free DOFs and generate no equation.
 *
 * @param model_data Nodal position field used to evaluate local coordinate axes.
 * @param equations Constraint collection receiving the generated rows.
 * @param node_id Global node identifier to constrain.
 */
void Support::apply_to_node(model::ModelData& model_data, constraint::Equations& equations, ID node_id) {
    // Validate the geometric field and node identifier before evaluating a
    // potentially position-dependent local coordinate system
    logging::error(model_data.positions != nullptr,
        "SUPPORT: positions field is not initialized");
    logging::error(node_id >= 0 && static_cast<Index>(node_id) < model_data.positions->rows,
        "SUPPORT: node ", node_id, " is outside the compiled node domain");

    // Evaluate the coordinate system at the global nodal position. The position
    // field may store generalized rows, but only the first three components are
    // geometric coordinates.
    const Vec6 position_vec = model_data.positions->row_vec6(static_cast<Index>(node_id));
    const Vec3 position     = position_vec.head<3>();

    // Traverse translations first and rotations second in the common structural
    // ordering [Ux, Uy, Uz, Rx, Ry, Rz]
    for (Dim dof = 0; dof < 6; ++dof) {
        // NaN is the semantic marker for an unconstrained component
        if (std::isnan(values_[dof])) {
            continue;
        }

        if (!coordinate_system_) {
            // A global prescription is one row of the identity mapping:
            //
            //     1 * u_dof = value_dof
            const constraint::EquationEntry entry{node_id, dof, Precision(1)};
            equations.emplace_back(
                std::initializer_list<constraint::EquationEntry>{entry},
                values_[dof]
            );
            continue;
        }

        // Evaluate the local orthonormal basis at this node. `get_axes()` uses
        // the coordinate system's own local parameterization, so the global
        // position is mapped first.
        const Vec3 local_position = coordinate_system_->to_local(position);
        const auto axes           = coordinate_system_->get_axes(local_position);

        // The local axis index is shared by translational and rotational blocks.
        // The block offset selects [Ux,Uy,Uz] or [Rx,Ry,Rz] in global numbering.
        const Dim  axis         = dof % 3;
        const Dim  block_offset = (dof / 3) * 3;
        const Vec3 direction    = axes.col(axis);

        // Expand e_a^T u = value into the three corresponding global DOFs
        const constraint::EquationEntry entry_x{node_id, static_cast<Dim>(block_offset + 0), direction[0]};
        const constraint::EquationEntry entry_y{node_id, static_cast<Dim>(block_offset + 1), direction[1]};
        const constraint::EquationEntry entry_z{node_id, static_cast<Dim>(block_offset + 2), direction[2]};

        equations.emplace_back(
            std::initializer_list<constraint::EquationEntry>{entry_x, entry_y, entry_z},
            values_[dof]
        );
    }
}

/**
 * Builds a compact diagnostic representation of the support definition.
 *
 * The output reports the semantic target, every finite prescribed generalized
 * component and the optional local coordinate system. The representation uses
 * the stored definition only and does not expand element or surface connectivity.
 *
 * @return Human-readable support description.
 */
std::string Support::str() const {
    // Resolve the one active semantic target to a compact type-qualified label
    const auto region_name = [&]() -> std::string {
        if (node_region_)
            return "NSET " + node_region_->name + " (" + std::to_string(node_region_->size()) + ")";
        if (surface_region_)
            return "SFSET " + surface_region_->name + " (" + std::to_string(surface_region_->size()) + ")";
        if (element_region_)
            return "ELSET " + element_region_->name + " (" + std::to_string(element_region_->size()) + ")";
        return "(unknown)";
    }();

    // Labels follow the generalized storage order in `values_`
    static constexpr const char* labels[] = {
        "Ux", "Uy", "Uz", "Rx", "Ry", "Rz"
    };

    std::ostringstream os;
    os << "Support: target=" << region_name << ", dof=[";

    // Print only finite prescriptions and avoid a trailing separator
    bool first = true;
    for (Dim dof = 0; dof < 6; ++dof) {
        if (std::isnan(values_[dof])) {
            continue;
        }

        if (!first) {
            os << ", ";
        }

        os << labels[dof] << '=' << values_[dof];
        first = false;
    }

    os << ']';

    if (coordinate_system_) {
        os << ", orientation=" << coordinate_system_->name;
    }

    return os.str();
}

} // namespace bc
} // namespace fem
