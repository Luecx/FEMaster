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

namespace fem {
namespace constraint {

/**
 * @copydoc Tie::Tie
 */
Tie::Tie(model::SurfaceRegion::Ptr master,
         model::NodeRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_set(std::move(master))
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
        Vec2 best_local;

        for (ID s_id : *master_set) {
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

        if (best_dist > distance) {
            continue;
        }

        auto s_ptr = surfaces.at(best_id);

        if (adjust) {
            auto mapped = s_ptr->local_to_global(best_local, node_coords);
            node_coords(id, 0) = mapped(0);
            node_coords(id, 1) = mapped(1);
            node_coords(id, 2) = mapped(2);
        }

        Dofs dofs_mask;
        for (int i = 0; i < 6; ++i) {
            dofs_mask(0, i) = system_nodal_dofs(id, i) >= 0;
        }
        for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
            ID master_node_id = s_ptr->nodes()[local_id];
            for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
                dofs_mask(0, dof_id) = dofs_mask(0, dof_id) && (system_nodal_dofs(master_node_id, dof_id) >= 0);
            }
        }

        auto nodal_contributions = s_ptr->shape_function(best_local);

        logging::warning(false,
                         nodal_contributions.array().abs().maxCoeff() < 10,
                         "Nodal contributions for slave node ",
                         id,
                         " exceed 10.");
        logging::error(false,
                       nodal_contributions.array().abs().maxCoeff() < 100,
                       "Nodal contributions for slave node ",
                       id,
                       " exceed 100.");

        for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
            if (!dofs_mask(0, dof_id)) {
                continue;
            }

            std::vector<EquationEntry> entries;
            entries.reserve(static_cast<std::size_t>(s_ptr->n_nodes) + 1);
            entries.emplace_back(EquationEntry{id, dof_id, Precision(1)});

            for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
                ID master_node_id = s_ptr->nodes()[local_id];
                Precision weight = nodal_contributions(local_id);
                entries.emplace_back(EquationEntry{master_node_id, dof_id, -weight});
            }

            equations.emplace_back(entries, Precision(0));
        }
    }

    return equations;
}

} // namespace constraint
} // namespace fem
