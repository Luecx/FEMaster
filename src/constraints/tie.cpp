/**
 * @file tie.cpp
 * @brief Implements the tie constraint used to couple slave nodes to master surfaces or lines.
 *
 * The implementation projects slave nodes onto candidate master geometries and
 * assembles the corresponding compatibility equations.
 *
 * Slave sets can be provided as:
 *  - node sets (direct node IDs)
 *  - surface sets (nodes are extracted from the referenced surfaces)
 *
 * @see src/constraints/tie.h
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "tie.h"

#include "../model/model_data.h"

#include <limits>
#include <unordered_set>
#include <vector>
#include "bvh.h"

namespace fem {
namespace constraint {

Tie::Tie(model::SurfaceRegion::Ptr master,
         model::NodeRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(std::move(master))
    , master_lines(nullptr)
    , slave_nodes(std::move(slave))
    , slave_surfaces(nullptr)
    , distance(max_distance)
    , adjust(do_adjust) {}

Tie::Tie(model::SurfaceRegion::Ptr master,
         model::SurfaceRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(std::move(master))
    , master_lines(nullptr)
    , slave_nodes(nullptr)
    , slave_surfaces(std::move(slave))
    , distance(max_distance)
    , adjust(do_adjust) {}

Tie::Tie(model::LineRegion::Ptr master,
         model::NodeRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(nullptr)
    , master_lines(std::move(master))
    , slave_nodes(std::move(slave))
    , slave_surfaces(nullptr)
    , distance(max_distance)
    , adjust(do_adjust) {}

Tie::Tie(model::LineRegion::Ptr master,
         model::SurfaceRegion::Ptr slave,
         Precision max_distance,
         bool do_adjust)
    : master_surfaces(nullptr)
    , master_lines(std::move(master))
    , slave_nodes(nullptr)
    , slave_surfaces(std::move(slave))
    , distance(max_distance)
    , adjust(do_adjust) {}

static std::vector<ID> collect_slave_nodes(const model::NodeRegion::Ptr& slave_nodes,
                                           const model::SurfaceRegion::Ptr& slave_surfaces,
                                           model::ModelData& model_data) {
    std::vector<ID> out;

    if (slave_nodes) {
        out.reserve(static_cast<std::size_t>(slave_nodes->size()));
        for (ID id : *slave_nodes) {
            out.push_back(id);
        }
        return out;
    }

    // Surface-slave: extract nodes from each surface, unique them.
    std::unordered_set<ID> unique;
    auto& surfaces = model_data.surfaces;

    for (ID s_id : *slave_surfaces) {
        if (static_cast<std::size_t>(s_id) >= surfaces.size()) {
            continue;
        }

        auto s_ptr = surfaces[static_cast<std::size_t>(s_id)];
        if (s_ptr == nullptr) {
            continue;
        }

        for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
            unique.insert(s_ptr->nodes()[local_id]);
        }
    }

    out.reserve(unique.size());
    for (ID nid : unique) {
        out.push_back(nid);
    }

    return out;
}

/**
 * @copydoc Tie::get_equations
 */
Equations Tie::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    auto& node_coords = *model_data.positions;
    auto& surfaces    = model_data.surfaces;
    auto& lines       = model_data.lines;

    Equations equations;

	// build the bvh for fast access to elements to not query all elements for every slave node
	BvhAabb bvh(distance);
	if (master_surfaces) {
        for (ID s_id : *master_surfaces) {
			if (static_cast<std::size_t>(s_id) >= surfaces.size()) {
                continue;
            }
            auto s_ptr = surfaces[static_cast<std::size_t>(s_id)];
            if (s_ptr == nullptr) {
                continue;
            }
			bvh.add_element(s_id, node_coords, s_ptr->nodes(), s_ptr->n_nodes);
        }
    } else if (master_lines) {
        for (ID l_id : *master_lines) {
            if (static_cast<std::size_t>(l_id) >= lines.size()) {
                continue;
            }
            auto l_ptr = lines[static_cast<std::size_t>(l_id)];
            if (l_ptr == nullptr) {
                continue;
            }
			bvh.add_element(l_id, node_coords, l_ptr->nodes(), l_ptr->n_nodes);
		}
	}
	bvh.finalize();

    // Build the actual list of slave node IDs (direct node-set or extracted from slave surface-set).
    std::vector<ID> slave_node_ids = collect_slave_nodes(slave_nodes, slave_surfaces, model_data);

	std::vector<ID> candidates;
	candidates.reserve(64);

    for (ID id : slave_node_ids) {
        // ---------------------------------------------------------------------
        // Slave node position
        // ---------------------------------------------------------------------
        const Index node_idx = static_cast<Index>(id);
        Vec3 node_pos = node_coords.row_vec3(node_idx);

        // ---------------------------------------------------------------------
        // Find closest master geometry (surface or line)
        // ---------------------------------------------------------------------
        Precision best_dist = std::numeric_limits<Precision>::max();
        ID        best_id   = -1;
        Vec2      best_local; // surfaces: (r,s); lines: (r,0)

        const auto& cand_ids = bvh.query_point(node_pos, &candidates);

        if (master_surfaces) {
            for (ID s_id : cand_ids) {
                if (static_cast<std::size_t>(s_id) >= surfaces.size()) {
                    continue;
                }

                auto s_ptr = surfaces[static_cast<std::size_t>(s_id)];
                if (s_ptr == nullptr) {
                    continue;
                }

                const Vec2 local     = s_ptr->global_to_local(node_pos, node_coords, true);
                const Vec3 mapped    = s_ptr->local_to_global(local, node_coords);
                const Precision dist = (node_pos - mapped).norm();

                if (dist > distance) {
                    continue;
                }

                if (dist < best_dist) {
                    best_id    = s_id;
                    best_dist  = dist;
                    best_local = local;
                }
            }
        } else if (master_lines) {
            for (ID l_id : cand_ids) {
                if (static_cast<std::size_t>(l_id) >= lines.size()) {
                    continue;
                }

                auto l_ptr = lines[static_cast<std::size_t>(l_id)];
                if (l_ptr == nullptr) {
                    continue;
                }

                const Precision r     = l_ptr->global_to_local(node_pos, node_coords, true);
                const Vec3 mapped     = l_ptr->local_to_global(r, node_coords);
                const Precision dist  = (node_pos - mapped).norm();

                if (dist > distance) {
                    continue;
                }

                if (dist < best_dist) {
                    best_id    = l_id;
                    best_dist  = dist;
                    best_local = Vec2(r, Precision(0));
                }
            }
        }

        if (best_id < 0 || best_dist > distance) {
            continue;
        }

        // ---------------------------------------------------------------------
        // Map back to global and optionally snap (adjust slave nodes)
        // ---------------------------------------------------------------------
        Vec3 mapped_pos = Vec3::Zero();
        if (master_surfaces) {
            auto s_ptr = surfaces[static_cast<std::size_t>(best_id)];
            mapped_pos = s_ptr->local_to_global(best_local, node_coords);
        } else {
            auto l_ptr = lines[static_cast<std::size_t>(best_id)];
            mapped_pos = l_ptr->local_to_global(best_local(0), node_coords);
        }

        if (adjust) {
            node_coords(node_idx, 0) = mapped_pos(0);
            node_coords(node_idx, 1) = mapped_pos(1);
            node_coords(node_idx, 2) = mapped_pos(2);
        }

        // ---------------------------------------------------------------------
        // Determine active DOFs: slave must have DOF AND all master nodes must have DOF
        // ---------------------------------------------------------------------
        Dofs dofs_mask;
        for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
            dofs_mask(0, dof_id) = (system_nodal_dofs(id, dof_id) >= 0);
        }

        DynamicVector nodal_contributions;
        Index n_master_nodes = 0;

        if (master_surfaces) {
            auto s_ptr = surfaces[static_cast<std::size_t>(best_id)];

            for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
                const ID master_node_id = s_ptr->nodes()[local_id];
                for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
                    dofs_mask(0, dof_id) =
                        dofs_mask(0, dof_id) && (system_nodal_dofs(master_node_id, dof_id) >= 0);
                }
            }

            nodal_contributions = s_ptr->shape_function(best_local);
            n_master_nodes      = s_ptr->n_nodes;
        } else {
            auto l_ptr = lines[static_cast<std::size_t>(best_id)];

            for (ID local_id = 0; local_id < static_cast<ID>(l_ptr->n_nodes); ++local_id) {
                const ID master_node_id = l_ptr->nodes()[local_id];
                for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
                    dofs_mask(0, dof_id) =
                        dofs_mask(0, dof_id) && (system_nodal_dofs(master_node_id, dof_id) >= 0);
                }
            }

            nodal_contributions = l_ptr->shape_function(best_local(0));
            n_master_nodes      = l_ptr->n_nodes;
        }

        logging::warning(
            nodal_contributions.array().abs().maxCoeff() < 10,
            "Nodal contributions for slave node ", id, " exceed 10.");
        logging::error(
            nodal_contributions.array().abs().maxCoeff() < 100,
            "Nodal contributions for slave node ", id, " exceed 100.");

        // ---------------------------------------------------------------------
        // Assemble compatibility equations: u_slave - sum_i N_i u_master_i = 0 per DOF
        // ---------------------------------------------------------------------
        for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
            if (!dofs_mask(0, dof_id)) {
                continue;
            }

            std::vector<EquationEntry> entries;
            entries.reserve(static_cast<std::size_t>(n_master_nodes) + 1);
            entries.emplace_back(EquationEntry{id, dof_id, Precision(1)});

            if (master_surfaces) {
                auto s_ptr = surfaces[static_cast<std::size_t>(best_id)];
                for (ID local_id = 0; local_id < static_cast<ID>(s_ptr->n_nodes); ++local_id) {
                    const ID master_node_id = s_ptr->nodes()[local_id];
                    const Precision w        = nodal_contributions(local_id);
                    entries.emplace_back(EquationEntry{master_node_id, dof_id, -w});
                }
            } else {
                auto l_ptr = lines[static_cast<std::size_t>(best_id)];
                for (ID local_id = 0; local_id < static_cast<ID>(l_ptr->n_nodes); ++local_id) {
                    const ID master_node_id = l_ptr->nodes()[local_id];
                    const Precision w        = nodal_contributions(local_id);
                    entries.emplace_back(EquationEntry{master_node_id, dof_id, -w});
                }
            }

            equations.emplace_back(entries, Precision(0));
        }
    }

    return equations;
}

} // namespace constraint
} // namespace fem
