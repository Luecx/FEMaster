/**
 * @file contact.cpp
 * @brief Implements simple frictionless node-to-surface penalty contact.
 */

#include "contact.h"

#include "../../model/model_data.h"
#include "bvh.h"

#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fem {
namespace constraint {
namespace {

struct DofCoeff {
    int       dof_id;
    Precision coeff;
};

std::vector<ID> collect_slave_nodes(const model::NodeRegion::Ptr& slave_nodes,
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

void add_translational_force(model::NodeData& nodal_forces,
                             ID node_id,
                             const Vec3& force) {
    const Index row = static_cast<Index>(node_id);
    for (Dim dof = 0; dof < 3; ++dof) {
        nodal_forces(row, dof) += force(dof);
    }
}

void add_gap_coefficients(std::vector<DofCoeff>& coeffs,
                          SystemDofIds& system_nodal_dofs,
                          ID node_id,
                          const Vec3& direction,
                          Precision weight) {
    for (Dim dof = 0; dof < 3; ++dof) {
        const int global_dof = system_nodal_dofs(node_id, dof);
        if (global_dof < 0) {
            continue;
        }

        const Precision coeff = weight * direction(dof);
        if (std::abs(coeff) <= Precision(1e-14)) {
            continue;
        }

        coeffs.push_back({global_dof, coeff});
    }
}

} // namespace

Contact::Contact(model::SurfaceRegion::Ptr master,
                 model::NodeRegion::Ptr slave,
                 Precision search_distance,
                 Precision penalty_stiffness,
                 Precision contact_clearance,
                 bool flip_master_normal)
    : master_surfaces(std::move(master))
    , slave_nodes(std::move(slave))
    , slave_surfaces(nullptr)
    , distance(search_distance)
    , penalty(penalty_stiffness)
    , clearance(contact_clearance)
    , flip_normal(flip_master_normal) {
    logging::error(distance >= Precision(0), "CONTACT: DISTANCE must be non-negative");
    logging::error(penalty > Precision(0), "CONTACT: PENALTY must be positive");
}

Contact::Contact(model::SurfaceRegion::Ptr master,
                 model::SurfaceRegion::Ptr slave,
                 Precision search_distance,
                 Precision penalty_stiffness,
                 Precision contact_clearance,
                 bool flip_master_normal)
    : master_surfaces(std::move(master))
    , slave_nodes(nullptr)
    , slave_surfaces(std::move(slave))
    , distance(search_distance)
    , penalty(penalty_stiffness)
    , clearance(contact_clearance)
    , flip_normal(flip_master_normal) {
    logging::error(distance >= Precision(0), "CONTACT: DISTANCE must be non-negative");
    logging::error(penalty > Precision(0), "CONTACT: PENALTY must be positive");
}

void Contact::assemble(SystemDofIds& system_nodal_dofs,
                       model::ModelData& model_data,
                       model::NodeData& nodal_forces,
                       TripletList& triplets) const {
    logging::error(model_data.positions != nullptr, "CONTACT: positions field not set in model data");

    auto& node_coords = *model_data.positions;
    auto& surfaces    = model_data.surfaces;

    BvhAabb bvh(distance);
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
    bvh.finalize();

    if (!bvh.valid()) {
        return;
    }

    std::vector<ID> slave_node_ids = collect_slave_nodes(slave_nodes, slave_surfaces, model_data);
    std::vector<ID> candidates;
    candidates.reserve(64);

    std::vector<DofCoeff> gap_coeffs;
    gap_coeffs.reserve(27);

    for (ID slave_node_id : slave_node_ids) {
        const Vec3 slave_pos = node_coords.row_vec3(static_cast<Index>(slave_node_id));

        Precision best_distance = std::numeric_limits<Precision>::max();
        ID        best_surface  = -1;
        Vec2      best_local    = Vec2::Zero();

        const auto& cand_ids = bvh.query_point(slave_pos, &candidates);
        for (ID s_id : cand_ids) {
            if (static_cast<std::size_t>(s_id) >= surfaces.size()) {
                continue;
            }

            auto s_ptr = surfaces[static_cast<std::size_t>(s_id)];
            if (s_ptr == nullptr) {
                continue;
            }

            const Vec2 local = s_ptr->global_to_local(slave_pos, node_coords, true);
            const Vec3 mapped = s_ptr->local_to_global(local, node_coords);
            const Precision distance_to_surface = (slave_pos - mapped).norm();

            if (distance_to_surface > distance) {
                continue;
            }

            if (distance_to_surface < best_distance) {
                best_distance = distance_to_surface;
                best_surface  = s_id;
                best_local    = local;
            }
        }

        if (best_surface < 0) {
            continue;
        }

        auto surface = surfaces[static_cast<std::size_t>(best_surface)];
        const Vec3 master_pos = surface->local_to_global(best_local, node_coords);
        Vec3 normal = surface->normal(node_coords, best_local);
        if (flip_normal) {
            normal = -normal;
        }

        logging::error(normal.allFinite(),
                       "CONTACT: surface ", best_surface, " has invalid normal");

        const Precision gap = (slave_pos - master_pos).dot(normal) - clearance;
        if (gap >= Precision(0)) {
            continue;
        }

        const DynamicVector shape = surface->shape_function(best_local);
        const Vec3 slave_force = (penalty * gap) * normal;

        add_translational_force(nodal_forces, slave_node_id, slave_force);

        for (ID local_id = 0; local_id < static_cast<ID>(surface->n_nodes); ++local_id) {
            const ID master_node_id = surface->nodes()[local_id];
            const Vec3 master_force = -shape(local_id) * slave_force;
            add_translational_force(nodal_forces, master_node_id, master_force);
        }

        gap_coeffs.clear();
        add_gap_coefficients(gap_coeffs,
                             system_nodal_dofs,
                             slave_node_id,
                             normal,
                             Precision(1));

        for (ID local_id = 0; local_id < static_cast<ID>(surface->n_nodes); ++local_id) {
            const ID master_node_id = surface->nodes()[local_id];
            add_gap_coefficients(gap_coeffs,
                                 system_nodal_dofs,
                                 master_node_id,
                                 normal,
                                 -shape(local_id));
        }

        for (const auto& row : gap_coeffs) {
            for (const auto& col : gap_coeffs) {
                triplets.emplace_back(row.dof_id,
                                      col.dof_id,
                                      penalty * row.coeff * col.coeff);
            }
        }
    }
}

} // namespace constraint
} // namespace fem
