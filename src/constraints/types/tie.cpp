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
 * @see src/constraints/types/tie.h
 * @see src/constraints/types/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "tie.h"

#include "../../model/model_data.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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

Equations Tie::get_surface_surface_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    logging::error(!adjust,
                   "TIE: ADJUST is not supported for surface-surface mortar");
    logging::error(model_data.positions != nullptr,
                   "TIE: positions field not set in model data");

    constexpr Precision tolerance = Precision(1e-12);

    auto& node_coords = *model_data.positions;
    auto& surfaces    = model_data.surfaces;

    struct MasterPatch {
        ID                 surface_id{};
        std::array<Vec2, 3> local{};
        std::array<Vec3, 3> global{};
        BvhAabb::Aabb      box{};
    };

    struct MortarRow {
        std::unordered_map<ID, Precision> slave_coeffs{};
        std::unordered_map<ID, Precision> master_coeffs{};
    };

    auto cross_2d = [](const Vec2& a, const Vec2& b) {
        return a(0) * b(1) - a(1) * b(0);
    };

    auto local_triangles = [](const model::SurfaceInterface::Ptr& surface) {
        std::vector<std::array<Vec2, 3>> triangles{};

        if (surface->n_nodes == 6) {
            const Vec2 p0(Precision(0.0), Precision(0.0));
            const Vec2 p1(Precision(1.0), Precision(0.0));
            const Vec2 p2(Precision(0.0), Precision(1.0));
            const Vec2 p3(Precision(0.5), Precision(0.0));
            const Vec2 p4(Precision(0.5), Precision(0.5));
            const Vec2 p5(Precision(0.0), Precision(0.5));

            return std::vector<std::array<Vec2, 3>>{
                {p0, p3, p5},
                {p3, p1, p4},
                {p5, p4, p2},
                {p3, p4, p5}
            };
        }

        if (surface->n_nodes == 8) {
            const Vec2 center(Precision(0.0), Precision(0.0));
            const std::array<Vec2, 8> boundary{
                Vec2(Precision(-1.0), Precision(-1.0)),
                Vec2(Precision( 0.0), Precision(-1.0)),
                Vec2(Precision( 1.0), Precision(-1.0)),
                Vec2(Precision( 1.0), Precision( 0.0)),
                Vec2(Precision( 1.0), Precision( 1.0)),
                Vec2(Precision( 0.0), Precision( 1.0)),
                Vec2(Precision(-1.0), Precision( 1.0)),
                Vec2(Precision(-1.0), Precision( 0.0))
            };

            triangles.reserve(boundary.size());
            for (std::size_t i = 0; i < boundary.size(); ++i) {
                triangles.push_back({boundary[i], boundary[(i + 1) % boundary.size()], center});
            }

            return triangles;
        }

        const auto domain = surface->local_domain_polygon();

        if (domain.size() < 3) {
            return triangles;
        }

        triangles.reserve(domain.size() - 2);
        for (std::size_t i = 1; i + 1 < domain.size(); ++i) {
            triangles.push_back({domain[0], domain[i], domain[i + 1]});
        }

        return triangles;
    };

    auto make_surface_aabb = [&](const model::SurfaceInterface::Ptr& surface) {
        BvhAabb::Aabb box = BvhAabb::Aabb::invalid();

        for (ID local_id = 0; local_id < static_cast<ID>(surface->n_nodes); ++local_id) {
            const ID node_id = surface->nodes()[local_id];
            box.expand_point(node_coords.row_vec3(static_cast<Index>(node_id)));
        }

        return box;
    };

    auto make_master_patch = [&](ID surface_id,
                                 const model::SurfaceInterface::Ptr& surface,
                                 const std::array<Vec2, 3>& local) {
        MasterPatch patch;
        patch.surface_id = surface_id;
        patch.local      = local;
        patch.box        = BvhAabb::Aabb::invalid();

        for (std::size_t i = 0; i < patch.local.size(); ++i) {
            patch.global[i] = surface->local_to_global(patch.local[i], node_coords);
            patch.box.expand_point(patch.global[i]);
        }

        return patch;
    };

    auto barycentric = [&](const Vec2& point, const std::array<Vec2, 3>& triangle) {
        const Vec2 edge_r = triangle[1] - triangle[0];
        const Vec2 edge_s = triangle[2] - triangle[0];
        const Vec2 offset = point       - triangle[0];

        const Precision denominator = cross_2d(edge_r, edge_s);
        const Precision lambda_1    = cross_2d(offset, edge_s) / denominator;
        const Precision lambda_2    = cross_2d(edge_r, offset) / denominator;

        Vec3 lambda;
        lambda << Precision(1) - lambda_1 - lambda_2,
                  lambda_1,
                  lambda_2;
        return lambda;
    };

    std::vector<MasterPatch> master_patches{};
    BvhAabb master_bvh;

    for (ID master_surface_id : *master_surfaces) {
        if (static_cast<std::size_t>(master_surface_id) >= surfaces.size()) {
            continue;
        }

        auto master = surfaces[static_cast<std::size_t>(master_surface_id)];
        if (master == nullptr) {
            continue;
        }

        for (const auto& local_triangle : local_triangles(master)) {
            const ID patch_id = static_cast<ID>(master_patches.size());
            master_patches.push_back(make_master_patch(master_surface_id, master, local_triangle));
            master_bvh.add_aabb(patch_id, master_patches.back().box);
        }
    }

    master_bvh.finalize();

    Equations equations;
    if (!master_bvh.valid()) {
        return equations;
    }

    const math::quadrature::Quadrature triangle_quadrature{
        math::quadrature::DOMAIN_ISO_TRI,
        math::quadrature::ORDER_QUARTIC
    };

    std::vector<ID> candidate_patch_ids{};
    candidate_patch_ids.reserve(64);

    std::unordered_map<ID, MortarRow> mortar_rows{};

    for (ID slave_surface_id : *slave_surfaces) {
        if (static_cast<std::size_t>(slave_surface_id) >= surfaces.size()) {
            continue;
        }

        auto slave = surfaces[static_cast<std::size_t>(slave_surface_id)];
        if (slave == nullptr) {
            continue;
        }

        const Index n_slave = slave->n_nodes;

        Precision slave_area = Precision(0);

        for (const auto& local_triangle : local_triangles(slave)) {
            slave->integrate_triangular(
                node_coords,
                model::SurfaceInterface::Polygon{
                    local_triangle[0],
                    local_triangle[1],
                    local_triangle[2]
                },
                triangle_quadrature,
                [&](const Vec2&, const Vec3&, Precision weight) {
                    slave_area += weight;
                });
        }

        logging::error(slave_area > tolerance,
                       "TIE: slave surface ",
                       slave_surface_id,
                       " has zero mortar area");

        BvhAabb::Aabb slave_box = make_surface_aabb(slave);
        slave_box.inflate(distance);

        const auto& candidates = master_bvh.query_aabb(slave_box, &candidate_patch_ids);

        std::vector<ID> candidate_patches{};
        candidate_patches.reserve(candidates.size());

        for (ID patch_id : candidates) {
            if (patch_id < 0 || static_cast<std::size_t>(patch_id) >= master_patches.size()) {
                continue;
            }

            const MasterPatch& patch = master_patches[static_cast<std::size_t>(patch_id)];
            auto master = surfaces[static_cast<std::size_t>(patch.surface_id)];
            if (master == nullptr) {
                continue;
            }

            candidate_patches.push_back(patch_id);
        }

        if (candidate_patches.empty()) {
            logging::warning(false,
                             "TIE: dual mortar slave surface ",
                             slave_surface_id,
                             " has no master candidates");
            continue;
        }

        std::vector<ID> local_master_node_ids{};
        std::unordered_map<ID, Index> local_master_index{};

        for (ID patch_id : candidate_patches) {
            const MasterPatch& patch = master_patches[static_cast<std::size_t>(patch_id)];
            auto master = surfaces[static_cast<std::size_t>(patch.surface_id)];

            for (Index col = 0; col < master->n_nodes; ++col) {
                const ID master_node = master->nodes()[col];
                if (local_master_index.find(master_node) != local_master_index.end()) {
                    continue;
                }

                const Index local_col = static_cast<Index>(local_master_node_ids.size());
                local_master_node_ids.push_back(master_node);
                local_master_index[master_node] = local_col;
            }
        }

        const Index n_master = static_cast<Index>(local_master_node_ids.size());

        DynamicMatrix G = DynamicMatrix::Zero(n_slave, n_slave);
        DynamicMatrix B = DynamicMatrix::Zero(n_slave, n_master);

        for (ID patch_id : candidate_patches) {
            const MasterPatch& patch = master_patches[static_cast<std::size_t>(patch_id)];
            auto master = surfaces[static_cast<std::size_t>(patch.surface_id)];

            std::array<Vec2, 3> slave_local{};
            for (std::size_t i = 0; i < slave_local.size(); ++i) {
                slave_local[i] = slave->global_to_local(patch.global[i], node_coords, false);
            }

            const Vec2 edge_r = slave_local[1] - slave_local[0];
            const Vec2 edge_s = slave_local[2] - slave_local[0];
            if (std::abs(cross_2d(edge_r, edge_s)) <= tolerance) {
                continue;
            }

            slave->integrate_triangular(
                node_coords,
                model::SurfaceInterface::Polygon{
                    slave_local[0],
                    slave_local[1],
                    slave_local[2]
                },
                triangle_quadrature,
                [&](const Vec2& local, const Vec3& global, Precision weight) {
                    const Vec3 lambda = barycentric(local, slave_local);
                    const Vec2 master_local =
                        lambda(0) * patch.local[0] +
                        lambda(1) * patch.local[1] +
                        lambda(2) * patch.local[2];

                    const Vec3 master_global = master->local_to_global(master_local, node_coords);
                    if ((master_global - global).norm() > distance + tolerance) {
                        return;
                    }

                    const DynamicVector Ns = slave->shape_function(local);
                    const DynamicVector Nm = master->shape_function(master_local);

                    for (Index row = 0; row < n_slave; ++row) {
                        for (Index col = 0; col < n_slave; ++col) {
                            G(row, col) += weight * Ns(row) * Ns(col);
                        }

                        for (Index col = 0; col < master->n_nodes; ++col) {
                            const ID master_node = master->nodes()[col];
                            const auto col_it    = local_master_index.find(master_node);

                            if (col_it == local_master_index.end()) {
                                continue;
                            }

                            B(row, col_it->second) += weight * Ns(row) * Nm(col);
                        }
                    }
                });
        }

        auto decomposition = G.ldlt();
        if (decomposition.info() != Eigen::Success) {
            logging::warning(false,
                             "TIE: failed to factorize local mortar matrix for slave surface ",
                             slave_surface_id);
            continue;
        }

        const DynamicMatrix W = decomposition.solve(B);
        if (!W.allFinite()) {
            logging::warning(false,
                             "TIE: invalid local mortar coefficients for slave surface ",
                             slave_surface_id);
            continue;
        }

        for (Index row = 0; row < n_slave; ++row) {
            const ID row_node = slave->nodes()[row];

            Precision row_sum = Precision(0);
            for (Index col = 0; col < n_master; ++col) {
                row_sum += W(row, col);
            }

            if (std::abs(row_sum) <= tolerance) {
                continue;
            }

            const Precision row_weight = G(row, row);
            if (std::abs(row_weight) <= tolerance) {
                continue;
            }

            mortar_rows[row_node].slave_coeffs[row_node] += row_weight;

            for (Index col = 0; col < n_master; ++col) {
                mortar_rows[row_node].master_coeffs[local_master_node_ids[col]] += row_weight * W(row, col);
            }
        }
    }

    constexpr Precision coefficient_tolerance = Precision(1.6e-2);

    for (const auto& [slave_node_id, row] : mortar_rows) {
        std::unordered_map<ID, Precision> coefficients{};
        coefficients.reserve(row.slave_coeffs.size() + row.master_coeffs.size());

        Precision slave_sum  = Precision(0);
        Precision master_sum = Precision(0);
        Precision scale      = Precision(0);

        for (const auto& [node_id, value] : row.slave_coeffs) {
            slave_sum += value;
            scale      = std::max(scale, std::abs(value));
        }

        for (const auto& [node_id, value] : row.master_coeffs) {
            (void) node_id;
            scale = std::max(scale, std::abs(value));
        }

        if (scale <= tolerance) {
            continue;
        }

        for (const auto& [node_id, value] : row.slave_coeffs) {
            coefficients[node_id] += value;
        }

        for (const auto& [node_id, value] : row.master_coeffs) {
            if (std::abs(value / scale) <= coefficient_tolerance) {
                continue;
            }

            master_sum += value;
        }

        if (std::abs(master_sum) <= tolerance) {
            continue;
        }

        const Precision master_correction = slave_sum / master_sum;

        for (const auto& [node_id, value] : row.master_coeffs) {
            if (std::abs(value / scale) <= coefficient_tolerance) {
                continue;
            }

            coefficients[node_id] -= master_correction * value;
        }

        Precision row_scale = Precision(0);
        for (const auto& [node_id, value] : coefficients) {
            (void) node_id;
            row_scale = std::max(row_scale, std::abs(value));
        }

        if (row_scale <= tolerance) {
            continue;
        }

        for (Dim dof_id = 0; dof_id < 6; ++dof_id) {
            if (system_nodal_dofs(slave_node_id, dof_id) < 0) {
                continue;
            }

            std::vector<EquationEntry> entries{};
            entries.reserve(coefficients.size());

            bool dof_supported = true;
            bool has_master    = false;

            for (const auto& [node_id, value] : coefficients) {
                const Precision scaled_value = value / row_scale;
                if (std::abs(value / scale) <= coefficient_tolerance) {
                    continue;
                }

                if (system_nodal_dofs(node_id, dof_id) < 0) {
                    dof_supported = false;
                    break;
                }

                if (node_id != slave_node_id) {
                    has_master = true;
                }

                entries.emplace_back(EquationEntry{node_id, dof_id, scaled_value});
            }

            if (!dof_supported || !has_master) {
                continue;
            }

            equations.emplace_back(entries, Precision(0));
        }
    }

    return equations;
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

    if (master_surfaces && slave_surfaces) {
        return get_surface_surface_equations(system_nodal_dofs, model_data);
    }

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
                    if (std::abs(w) > 1e-12)
                        entries.emplace_back(EquationEntry{master_node_id, dof_id, -w});
                }
            } else {
                auto l_ptr = lines[static_cast<std::size_t>(best_id)];
                for (ID local_id = 0; local_id < static_cast<ID>(l_ptr->n_nodes); ++local_id) {
                    const ID master_node_id = l_ptr->nodes()[local_id];
                    const Precision w        = nodal_contributions(local_id);
                    if (std::abs(w) > 1e-12)
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
