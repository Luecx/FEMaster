/**
 * @file tie.cpp
 * @brief Implements the tie constraint used to couple slave nodes to surfaces.
 *
 * The implementation projects slave nodes onto candidate master surfaces and
 * assembles the corresponding compatibility equations.
 *
 * @see src/constraints/tie.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "tie.h"

#include "../model/model_data.h"
#include "../model/geometry/line/line2a.h"

namespace fem {
namespace constraint {

/**
 * @copydoc Tie::Tie
 */
Tie::Tie(model::SurfaceRegion::Ptr master,
         model::NodeRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(std::move(master))
    , master_lines(nullptr)
    , slave_set(std::move(slave))
    , distance(max_distance)
    , adjust(do_adjust) {}

Tie::Tie(model::LineRegion::Ptr master,
         model::NodeRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(nullptr)
    , master_lines(std::move(master))
    , slave_set(std::move(slave))
    , distance(max_distance)
    , adjust(do_adjust) {}

/**
 * @copydoc Tie::get_equations
 */
Equations Tie::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    auto& node_coords = model_data.get(model::NodeDataEntries::POSITION);
    auto& surfaces = model_data.surfaces;

    Equations equations;

    for (ID id : *slave_set) {
        Vec3 node_pos;
        node_pos(0) = node_coords(id, 0);
        node_pos(1) = node_coords(id, 1);
        node_pos(2) = node_coords(id, 2);

        Precision best_dist = static_cast<Precision>(1e36);
        ID best_id = -1;
        Vec2 best_local; // For surfaces: (r,s); For lines: (r,0)

        // Branch on geometry kind: surfaces (2D) vs lines (1D)
        if (master_surfaces) {
            for (ID s_id : *master_surfaces) {
                auto s_ptr = surfaces.at(s_id);
                if (s_ptr == nullptr) {
                    continue;
                }

                auto local = s_ptr->global_to_local(node_pos, node_coords, true);
                auto mapped = s_ptr->local_to_global(local, node_coords);
                auto dist = (node_pos - mapped).norm();
                if (dist > distance) {
                    continue;
                }

                if (dist < best_dist) {
                    best_id = s_id;
                    best_dist = dist;
                    best_local = local;
                }
            }
        } else if (master_lines) {
            // Search across line masters using Line2A (2-node) representation
            for (ID l_id : *master_lines) {
                if (static_cast<std::size_t>(l_id) >= model_data.lines.size()) {
                    continue;
                }
                const auto nodes = model_data.lines[static_cast<std::size_t>(l_id)];
                fem::model::Line2A line(nodes);

                Precision r = line.global_to_local(node_pos, node_coords, true);
                Vec3 mapped = line.local_to_global(r, node_coords);
                auto dist = (node_pos - mapped).norm();
                if (dist > distance) {
                    continue;
                }

                if (dist < best_dist) {
                    best_id = l_id;
                    best_dist = dist;
                    best_local = Vec2(r, 0);
                }
            }
        }

        if (best_dist > distance) {
            continue;
        }

        // Map back to global and optionally snap
        Vec3 mapped_pos;
        if (master_surfaces) {
            auto s_ptr = surfaces.at(best_id);
            mapped_pos = s_ptr->local_to_global(best_local, node_coords);
        } else {
            const auto nodes = model_data.lines[static_cast<std::size_t>(best_id)];
            fem::model::Line2A line(nodes);
            mapped_pos = line.local_to_global(best_local(0), node_coords);
        }

        if (adjust) {
            node_coords(id, 0) = mapped_pos(0);
            node_coords(id, 1) = mapped_pos(1);
            node_coords(id, 2) = mapped_pos(2);
        }

        Dofs dofs_mask;
        for (int i = 0; i < 6; ++i) {
            dofs_mask(0, i) = system_nodal_dofs(id, i) >= 0;
        }
        DynamicVector nodal_contributions;
        Index n_master_nodes = 0;
        if (master_surfaces) {
            auto s_ptr = surfaces.at(best_id);
            for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
                ID master_node_id = s_ptr->nodes()[local_id];
                for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
                    dofs_mask(0, dof_id) = dofs_mask(0, dof_id) && (system_nodal_dofs(master_node_id, dof_id) >= 0);
                }
            }
            nodal_contributions = s_ptr->shape_function(best_local);
            n_master_nodes = s_ptr->n_nodes;
        } else {
            const auto nodes = model_data.lines[static_cast<std::size_t>(best_id)];
            fem::model::Line2A line(nodes);
            // Master DOF availability: both line nodes must have the DOF
            for (ID local_id = 0; local_id < 2; ++local_id) {
                ID master_node_id = nodes[static_cast<std::size_t>(local_id)];
                for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
                    dofs_mask(0, dof_id) = dofs_mask(0, dof_id) && (system_nodal_dofs(master_node_id, dof_id) >= 0);
                }
            }
            auto N = line.shape_function(best_local(0));
            nodal_contributions.resize(2);
            nodal_contributions(0) = N(0);
            nodal_contributions(1) = N(1);
            n_master_nodes = 2;
        }

        logging::warning(
            nodal_contributions.array().abs().maxCoeff() < 10,
            "Nodal contributions for slave node ", id, " exceed 10.");
        logging::error(
            nodal_contributions.array().abs().maxCoeff() < 100,
            "Nodal contributions for slave node ", id, " exceed 100.");

        for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
            if (!dofs_mask(0, dof_id)) {
                continue;
            }

            std::vector<EquationEntry> entries;
            entries.reserve(static_cast<std::size_t>(n_master_nodes) + 1);
            entries.emplace_back(EquationEntry{id, dof_id, Precision(1)});

            if (master_surfaces) {
                auto s_ptr = surfaces.at(best_id);
                for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
                    ID master_node_id = s_ptr->nodes()[local_id];
                    Precision weight = nodal_contributions(local_id);
                    entries.emplace_back(EquationEntry{master_node_id, dof_id, -weight});
                }
            } else {
                const auto nodes = model_data.lines[static_cast<std::size_t>(best_id)];
                for (ID local_id = 0; local_id < 2; ++local_id) {
                    ID master_node_id = nodes[static_cast<std::size_t>(local_id)];
                    Precision weight = nodal_contributions(local_id);
                    entries.emplace_back(EquationEntry{master_node_id, dof_id, -weight});
                }
            }

            equations.emplace_back(entries, Precision(0));
        }
    }

    return equations;
}

} // namespace constraint
} // namespace fem
