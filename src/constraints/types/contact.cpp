/**
 * @file contact.cpp
 * @brief Implements frictionless faceted node-to-surface penalty contact.
 *
 * The master partner is initialized by a bounded closest-point search and then
 * tracked locally across the connected master patch when a previous partner
 * exists. Contact projection, gap, force direction, and tangent consistently
 * use the geometrical normal of the selected master facet. No independent
 * patch-normal averaging is applied.
 *
 * The compression-only penalty law is piecewise linear. Open contact therefore
 * contributes exactly zero force and zero tangent. Slave surface sets use fixed
 * positive lumped tributary areas initialized during the first assembly.
 * Explicit slave node sets retain the historic FEMaster interpretation in
 * which PENALTY is a nodal stiffness.
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"
#include "bvh.h"

#include <Eigen/Cholesky>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
namespace fem::constraint {
namespace {

constexpr Precision maximum_projection_derivative_norm = Precision(5000);
constexpr Precision partner_switch_distance_tolerance = Precision(2.5e-5);

enum class ProjectionMode {
    Interior,
    Edge,
    Corner
};

struct ProjectionInfo {
    ProjectionMode mode           = ProjectionMode::Interior;
    Vec2           edge_direction = Vec2::Zero();
};

struct ContactDiagnostics {
    Index slave_nodes             = 0;
    Index candidate_surfaces      = 0;
    Index maximum_candidates      = 0;
    Index zero_candidate_slaves   = 0;
    Index valid_projections       = 0;
    Index invalid_projections     = 0;
    Index distance_rejections     = 0;
    Index no_partner              = 0;
    Index open_closest_partner    = 0;
    Index active_contacts         = 0;
    Index interior_consistent     = 0;
    Index interior_fallback       = 0;
    Index edge_consistent         = 0;
    Index edge_fallback           = 0;
    Index corner_direct           = 0;
    Index activations             = 0;
    Index deactivations           = 0;
    Index partner_switches        = 0;

    Precision maximum_penetration            = Precision(0);
    Precision maximum_closest_distance       = Precision(0);
    Precision maximum_slave_force            = Precision(0);
    Precision slave_force_squared_sum        = Precision(0);
    Precision minimum_surface_jacobian       = std::numeric_limits<Precision>::infinity();
    Precision minimum_interior_hessian_ratio = std::numeric_limits<Precision>::infinity();
    Precision minimum_edge_hessian_ratio     = std::numeric_limits<Precision>::infinity();
    Precision maximum_projection_derivative  = Precision(0);
    Precision maximum_raw_asymmetry          = Precision(0);
    Precision maximum_local_tangent_norm     = Precision(0);

    ID             worst_slave   = ID(-1);
    ID             worst_surface = ID(-1);
    Vec2           worst_local   = Vec2::Zero();
    ProjectionMode worst_mode    = ProjectionMode::Interior;
    Precision      worst_gap     = Precision(0);
    Precision      worst_distance = Precision(0);

    std::uint64_t signature = 1469598103934665603ull;
};

const char* projection_mode_name(ProjectionMode mode) {
    switch (mode) {
        case ProjectionMode::Interior:
            return "INTERIOR";
        case ProjectionMode::Edge:
            return "EDGE";
        case ProjectionMode::Corner:
            return "CORNER";
    }

    return "UNKNOWN";
}

bool projection_is_inside_closed_domain(
    Index       number_of_nodes,
    const Vec2& local
);

struct EdgeKey {
    ID first  = ID(-1);
    ID second = ID(-1);
    ID middle = ID(-1);

    bool operator==(const EdgeKey& other) const {
        return first  == other.first &&
               second == other.second &&
               middle == other.middle;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& key) const {
        std::size_t seed = std::hash<ID>{}(key.first);

        seed ^= std::hash<ID>{}(key.second) +
                std::size_t(0x9e3779b9) +
                (seed << 6) +
                (seed >> 2);

        seed ^= std::hash<ID>{}(key.middle) +
                std::size_t(0x9e3779b9) +
                (seed << 6) +
                (seed >> 2);

        return seed;
    }
};

struct PatchProjection {
    bool      valid      = false;
    ID        surface_id = ID(-1);
    Vec2      local      = Vec2::Zero();
    Precision distance   = std::numeric_limits<Precision>::max();
};

std::array<Index, 3> edge_local_nodes(
    const model::SurfaceInterface& surface,
    Index                          edge
) {
    if (surface.n_nodes == 3) {
        switch (edge) {
            case 0: return {0, 1, Index(-1)};
            case 1: return {1, 2, Index(-1)};
            default: return {2, 0, Index(-1)};
        }
    }

    if (surface.n_nodes == 6) {
        switch (edge) {
            case 0: return {0, 1, 3};
            case 1: return {1, 2, 4};
            default: return {2, 0, 5};
        }
    }

    if (surface.n_nodes == 4) {
        switch (edge) {
            case 0: return {0, 1, Index(-1)};
            case 1: return {1, 2, Index(-1)};
            case 2: return {2, 3, Index(-1)};
            default: return {3, 0, Index(-1)};
        }
    }

    switch (edge) {
        case 0: return {0, 1, 4};
        case 1: return {1, 2, 5};
        case 2: return {2, 3, 6};
        default: return {3, 0, 7};
    }
}

EdgeKey make_edge_key(
    const model::SurfaceInterface& surface,
    Index                          edge
) {
    const auto local_nodes =
        edge_local_nodes(surface, edge);

    ID first =
        surface.nodes()[local_nodes[0]];

    ID second =
        surface.nodes()[local_nodes[1]];

    const ID middle =
        local_nodes[2] >= 0
            ? surface.nodes()[local_nodes[2]]
            : ID(-1);

    if (second < first) {
        std::swap(first, second);
    }

    return {first, second, middle};
}

ID crossed_edge_index(
    const model::SurfaceInterface& surface,
    const Vec2&                    local
) {
    const Precision r = local(0);
    const Precision s = local(1);

    Precision largest_violation = Precision(1e-10);
    ID        edge              = ID(-1);

    auto update = [&](Precision violation, ID candidate_edge) {
        if (violation > largest_violation) {
            largest_violation = violation;
            edge              = candidate_edge;
        }
    };

    if (surface.n_nodes == 3 || surface.n_nodes == 6) {
        update(-s, 0);
        update(r + s - Precision(1), 1);
        update(-r, 2);
    } else {
        update(-Precision(1) - s, 0);
        update( r - Precision(1), 1);
        update( s - Precision(1), 2);
        update(-Precision(1) - r, 3);
    }

    return edge;
}

std::vector<std::array<ID, 4>> build_master_patch_topology(
    const model::SurfaceRegion::Ptr&                 master_surfaces,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    std::vector<std::array<ID, 4>> edge_neighbors(surfaces.size());

    for (auto& neighbors : edge_neighbors) {
        neighbors.fill(ID(-1));
    }

    std::unordered_map<EdgeKey, std::pair<ID, Index>, EdgeKeyHash> edge_owner;

    for (ID surface_id : *master_surfaces) {
        if (surface_id < 0 ||
            static_cast<std::size_t>(surface_id) >= surfaces.size()) {
            continue;
        }

        const auto& surface =
            surfaces[static_cast<std::size_t>(surface_id)];

        if (!surface) {
            continue;
        }

        for (Index edge = 0; edge < surface->n_edges; ++edge) {
            const EdgeKey key =
                make_edge_key(*surface, edge);

            const auto owner =
                edge_owner.find(key);

            if (owner == edge_owner.end()) {
                edge_owner.emplace(key, std::make_pair(surface_id, edge));
                continue;
            }

            const ID    neighbor_id   = owner->second.first;
            const Index neighbor_edge = owner->second.second;

            logging::error(
                edge_neighbors[static_cast<std::size_t>(neighbor_id)]
                                               [static_cast<std::size_t>(neighbor_edge)] < 0,
                "CONTACT: non-manifold master edge shared by more than two surfaces"
            );

            edge_neighbors[static_cast<std::size_t>(surface_id)]
                                           [static_cast<std::size_t>(edge)] =
                neighbor_id;

            edge_neighbors[static_cast<std::size_t>(neighbor_id)]
                                           [static_cast<std::size_t>(neighbor_edge)] =
                surface_id;
        }
    }

    return edge_neighbors;
}

bool valid_surface_id(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           surfaces[static_cast<std::size_t>(surface_id)] != nullptr;
}

PatchProjection project_on_surface(
    ID                                                surface_id,
    const Vec3&                                       slave_position,
    const model::Field&                               node_coords,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    PatchProjection projection;

    if (!valid_surface_id(surface_id, surfaces)) {
        return projection;
    }

    const auto& surface =
        surfaces[static_cast<std::size_t>(surface_id)];

    const Vec2 local =
        surface->global_to_local(
            slave_position,
            node_coords,
            true
        );

    if (!local.allFinite() ||
        !projection_is_inside_closed_domain(surface->n_nodes, local)) {
        return projection;
    }

    const Vec3 master_position =
        surface->local_to_global(local, node_coords);

    const Precision distance_to_surface =
        (slave_position - master_position).norm();

    if (!std::isfinite(distance_to_surface)) {
        return projection;
    }

    projection.valid      = true;
    projection.surface_id = surface_id;
    projection.local      = local;
    projection.distance   = distance_to_surface;

    return projection;
}

PatchProjection walk_master_patch(
    ID                                                start_surface_id,
    const Vec3&                                       slave_position,
    const model::Field&                               node_coords,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const std::vector<std::array<ID, 4>>&             edge_neighbors
) {
    ID current_surface_id = start_surface_id;

    std::unordered_set<ID> visited;
    visited.reserve(
        std::min<std::size_t>(
            edge_neighbors.size(),
            std::size_t(32)
        )
    );

    const Index maximum_steps =
        static_cast<Index>(
            std::min<std::size_t>(
                edge_neighbors.size(),
                std::size_t(64)
            )
        );

    for (Index step = 0; step < maximum_steps; ++step) {
        if (!valid_surface_id(current_surface_id, surfaces) ||
            visited.find(current_surface_id) != visited.end()) {
            break;
        }

        visited.insert(current_surface_id);

        const auto& surface =
            surfaces[static_cast<std::size_t>(current_surface_id)];

        const Vec2 unclipped_local =
            surface->global_to_local(
                slave_position,
                node_coords,
                false
            );

        if (unclipped_local.allFinite() &&
            projection_is_inside_closed_domain(surface->n_nodes, unclipped_local)) {
            return project_on_surface(
                current_surface_id,
                slave_position,
                node_coords,
                surfaces
            );
        }

        const ID crossed_edge =
            unclipped_local.allFinite()
                ? crossed_edge_index(*surface, unclipped_local)
                : ID(-1);

        if (crossed_edge < 0 ||
            static_cast<std::size_t>(current_surface_id) >=
                edge_neighbors.size()) {
            break;
        }

        const ID neighbor_surface_id =
            edge_neighbors[static_cast<std::size_t>(current_surface_id)]
                                   [static_cast<std::size_t>(crossed_edge)];

        if (!valid_surface_id(neighbor_surface_id, surfaces) ||
            visited.find(neighbor_surface_id) != visited.end()) {
            break;
        }

        current_surface_id = neighbor_surface_id;
    }

    // A failed walk deliberately returns no boundary fallback. The caller then
    // performs a bounded BVH search and compares all admissible candidates.
    return {};
}

void hash_value(
    std::uint64_t& signature,
    std::uint64_t  value
) {
    constexpr std::uint64_t prime = 1099511628211ull;

    for (Index byte = 0; byte < 8; ++byte) {
        signature ^= value & 0xffull;
        signature *= prime;
        value >>= 8;
    }
}

std::vector<ID> collect_slave_nodes(
    const model::NodeRegion::Ptr&    slave_nodes,
    const model::SurfaceRegion::Ptr& slave_surfaces,
    model::ModelData&                model_data
) {
    std::vector<ID> out;

    if (slave_nodes) {
        out.reserve(static_cast<std::size_t>(slave_nodes->size()));

        for (ID node_id : *slave_nodes) {
            out.push_back(node_id);
        }

        std::sort(out.begin(), out.end());
        return out;
    }

    std::unordered_set<ID> unique_nodes;
    auto&                  surfaces = model_data.surfaces;

    for (ID surface_id : *slave_surfaces) {
        if (surface_id < 0 ||
            static_cast<std::size_t>(surface_id) >= surfaces.size()) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

        if (!surface) {
            continue;
        }

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            unique_nodes.insert(surface->nodes()[local_node]);
        }
    }

    out.reserve(unique_nodes.size());

    for (ID node_id : unique_nodes) {
        out.push_back(node_id);
    }

    std::sort(out.begin(), out.end());
    return out;
}

std::unordered_map<ID, Precision> build_slave_tributary_areas(
    const model::SurfaceRegion::Ptr& slave_surfaces,
    model::ModelData&                model_data,
    const model::Field&              node_coords
) {
    std::unordered_map<ID, Precision> areas;

    if (!slave_surfaces) {
        return areas;
    }

    auto& surfaces = model_data.surfaces;

    for (ID surface_id : *slave_surfaces) {
        if (surface_id < 0 ||
            static_cast<std::size_t>(surface_id) >= surfaces.size()) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

        if (!surface || surface->n_nodes <= 0) {
            continue;
        }

        const Precision surface_area = surface->area(node_coords);

        logging::error(
            std::isfinite(surface_area) && surface_area > Precision(0),
            "CONTACT: invalid slave surface area for surface ",
            surface_id
        );

        // Positive lumping avoids negative tributary weights for quadratic
        // serendipity surfaces.
        const Precision nodal_area =
            surface_area / static_cast<Precision>(surface->n_nodes);

        for (Index local_node = 0;
             local_node < surface->n_nodes;
             ++local_node) {
            areas[surface->nodes()[local_node]] += nodal_area;
        }
    }

    return areas;
}

struct PenaltyLaw {
    Precision value      = Precision(0);
    Precision derivative = Precision(0);
};

PenaltyLaw evaluate_penalty_law(Precision gap) {
    return gap < Precision(0)
        ? PenaltyLaw{gap, Precision(1)}
        : PenaltyLaw{};
}

void add_translational_force(
    model::NodeData& nodal_forces,
    ID               node_id,
    const Vec3&      force
) {
    const Index row = static_cast<Index>(node_id);

    for (Dim component = 0; component < 3; ++component) {
        nodal_forces(row, component) += force(component);
    }
}

Mat3 skew(const Vec3& vector) {
    Mat3 matrix;

    matrix << Precision(0), -vector(2),     vector(1),
                 vector(2), Precision(0), -vector(0),
                -vector(1),     vector(0), Precision(0);

    return matrix;
}

bool projection_is_inside_closed_domain(
    Index       number_of_nodes,
    const Vec2& local
) {
    constexpr Precision tolerance = Precision(1e-9);

    if (number_of_nodes == 3 || number_of_nodes == 6) {
        const Precision r = local(0);
        const Precision s = local(1);

        return r >= -tolerance &&
               s >= -tolerance &&
               r + s <= Precision(1) + tolerance;
    }

    if (number_of_nodes == 4 || number_of_nodes == 8) {
        const Precision r = local(0);
        const Precision s = local(1);

        return r >= Precision(-1) - tolerance &&
               r <= Precision( 1) + tolerance &&
               s >= Precision(-1) - tolerance &&
               s <= Precision( 1) + tolerance;
    }

    return false;
}

ProjectionInfo classify_projection(
    Index       number_of_nodes,
    const Vec2& local
) {
    constexpr Precision edge_tolerance = Precision(1e-7);

    ProjectionInfo info;

    if (number_of_nodes == 3 || number_of_nodes == 6) {
        const Precision r = local(0);
        const Precision s = local(1);
        const Precision t = Precision(1) - r - s;

        const bool on_r_edge = std::abs(r) <= edge_tolerance;
        const bool on_s_edge = std::abs(s) <= edge_tolerance;
        const bool on_t_edge = std::abs(t) <= edge_tolerance;

        /*
         * A clipped projection on two natural boundaries is a true patch
         * corner. Treating it as an edge gives the Newton tangent a spurious
         * one-dimensional sliding direction along an edge endpoint and can
         * drive the line search into a non-smooth contact branch.
         */
        if ((on_r_edge ? 1 : 0) +
            (on_s_edge ? 1 : 0) +
            (on_t_edge ? 1 : 0) >= 2) {
            info.mode = ProjectionMode::Corner;
            return info;
        }

        if (on_r_edge) {
            info.mode           = ProjectionMode::Edge;
            info.edge_direction = Vec2(Precision(0), Precision(1));
            return info;
        }

        if (on_s_edge) {
            info.mode           = ProjectionMode::Edge;
            info.edge_direction = Vec2(Precision(1), Precision(0));
            return info;
        }

        if (on_t_edge) {
            info.mode           = ProjectionMode::Edge;
            info.edge_direction = Vec2(Precision(1), Precision(-1));
            return info;
        }

        info.mode = ProjectionMode::Interior;
        return info;
    }

    if (number_of_nodes == 4 || number_of_nodes == 8) {
        const Precision r = local(0);
        const Precision s = local(1);

        const bool on_r_min_edge = std::abs(r + Precision(1)) <= edge_tolerance;
        const bool on_r_max_edge = std::abs(r - Precision(1)) <= edge_tolerance;
        const bool on_s_min_edge = std::abs(s + Precision(1)) <= edge_tolerance;
        const bool on_s_max_edge = std::abs(s - Precision(1)) <= edge_tolerance;

        /*
         * A point on two quadrilateral boundaries is constrained to a corner.
         * The corner branch intentionally uses a zero projection derivative and
         * the patch-averaged corner normal, which is more conservative than an
         * artificial edge tangent through the endpoint.
         */
        const bool on_r_edge = on_r_min_edge || on_r_max_edge;
        const bool on_s_edge = on_s_min_edge || on_s_max_edge;

        if (on_r_edge && on_s_edge) {
            info.mode = ProjectionMode::Corner;
            return info;
        }

        if (on_r_min_edge || on_r_max_edge) {
            info.mode           = ProjectionMode::Edge;
            info.edge_direction = Vec2(Precision(0), Precision(1));
            return info;
        }

        if (on_s_min_edge || on_s_max_edge) {
            info.mode           = ProjectionMode::Edge;
            info.edge_direction = Vec2(Precision(1), Precision(0));
            return info;
        }

        info.mode = ProjectionMode::Interior;
        return info;
    }

    info.mode = ProjectionMode::Corner;
    return info;
}

template<Index N>
void scatter_contact_tangent(
    const std::array<ID, N + 1>&                    node_ids,
    const StaticMatrix<3 * (N + 1), 3 * (N + 1)>& tangent,
    SystemDofIds&                                   system_nodal_dofs,
    TripletList&                                    triplets
) {
    constexpr Index local_nodes = N + 1;
    constexpr Index local_dofs  = 3 * local_nodes;

    std::array<int, local_dofs> global_dofs{};

    for (Index local_node = 0; local_node < local_nodes; ++local_node) {
        const ID node_id = node_ids[static_cast<std::size_t>(local_node)];

        for (Dim component = 0; component < 3; ++component) {
            const Index local_dof = 3 * local_node + component;

            global_dofs[static_cast<std::size_t>(local_dof)] =
                system_nodal_dofs(node_id, component);
        }
    }

    for (Index local_row = 0; local_row < local_dofs; ++local_row) {
        const int global_row =
            global_dofs[static_cast<std::size_t>(local_row)];

        if (global_row < 0) {
            continue;
        }

        for (Index local_col = 0; local_col < local_dofs; ++local_col) {
            const int global_col =
                global_dofs[static_cast<std::size_t>(local_col)];

            if (global_col < 0) {
                continue;
            }

            const Precision value = tangent(local_row, local_col);

            if (std::abs(value) <= Precision(1e-14)) {
                continue;
            }

            triplets.emplace_back(global_row, global_col, value);
        }
    }
}

template<Index N>
void assemble_contact_pair(
    const model::Surface<N>& surface,
    ID                       slave_node_id,
    const Vec2&              master_local,
    const ProjectionInfo&    projection_info,
    Precision                penalty,
    Precision                clearance,
    Precision                slave_weight,
    bool                     flip_normal,
    SystemDofIds&            system_nodal_dofs,
    const model::Field&      node_coords,
    model::NodeData&         nodal_forces,
    TripletList&             triplets,
    ContactDiagnostics&      diagnostics
) {
    constexpr Index local_nodes = N + 1;
    constexpr Index local_dofs  = 3 * local_nodes;

    using LocalMatrix   = StaticMatrix<local_dofs, local_dofs>;
    using LocalRow      = StaticMatrix<1, local_dofs>;
    using LocalMap      = StaticMatrix<3, local_dofs>;
    using ProjectionMap = StaticMatrix<2, local_dofs>;

    const Precision r = master_local(0);
    const Precision s = master_local(1);

    const StaticMatrix<N, 1> shape =
        surface.shape_function(r, s);

    const StaticMatrix<N, 2> shape_derivative =
        surface.shape_derivative(r, s);

    const StaticMatrix<N, 3> shape_second_derivative =
        surface.shape_second_derivative(r, s);

    const StaticMatrix<N, 3> master_coordinates =
        surface.node_coords_global(node_coords);

    const Vec3 slave_position =
        node_coords.row_vec3(static_cast<Index>(slave_node_id));

    const Vec3 master_position =
        master_coordinates.transpose() * shape;

    const Vec3 difference =
        slave_position - master_position;

    const StaticMatrix<3, 2> first =
        master_coordinates.transpose() * shape_derivative;

    const Vec3 tangent_r = first.col(0);
    const Vec3 tangent_s = first.col(1);

    const StaticMatrix<3, 3> second =
        master_coordinates.transpose() * shape_second_derivative;

    const Vec3 tangent_rr = second.col(0);
    const Vec3 tangent_ss = second.col(1);
    const Vec3 tangent_rs = second.col(2);

    const Vec3 normal_unnormalized =
        tangent_r.cross(tangent_s);

    const Precision surface_jacobian =
        normal_unnormalized.norm();

    diagnostics.minimum_surface_jacobian =
        std::min(
            diagnostics.minimum_surface_jacobian,
            surface_jacobian
        );

    logging::error(
        std::isfinite(surface_jacobian) &&
        surface_jacobian > Precision(1e-14),
        "CONTACT: singular master-surface Jacobian"
    );

    const Vec3 normal_base =
        normal_unnormalized / surface_jacobian;

    const Vec3 normal =
        flip_normal ? -normal_base : normal_base;

    const Precision gap =
        difference.dot(normal) - clearance;

    std::array<ID, local_nodes> local_node_ids{};
    local_node_ids[0] = slave_node_id;

    for (Index master_node = 0; master_node < N; ++master_node) {
        local_node_ids[static_cast<std::size_t>(master_node + 1)] =
            surface.nodes()[master_node];
    }

    const PenaltyLaw contact_law =
        evaluate_penalty_law(gap);

    if (gap >= Precision(0)) {
        return;
    }

    const Precision contact_scale =
        slave_weight * penalty;

    const Vec3 slave_force =
        contact_scale * contact_law.value * normal;

    add_translational_force(
        nodal_forces,
        slave_node_id,
        slave_force
    );

    for (Index master_node = 0; master_node < N; ++master_node) {
        add_translational_force(
            nodal_forces,
            surface.nodes()[master_node],
            -shape(master_node) * slave_force
        );
    }

    LocalMap direct_difference = LocalMap::Zero();
    direct_difference.template block<3, 3>(0, 0) = Mat3::Identity();

    LocalMap direct_tangent_r = LocalMap::Zero();
    LocalMap direct_tangent_s = LocalMap::Zero();

    for (Index master_node = 0; master_node < N; ++master_node) {
        const Index column_offset = 3 * (master_node + 1);

        direct_difference.template block<3, 3>(0, column_offset) =
            -shape(master_node) * Mat3::Identity();

        direct_tangent_r.template block<3, 3>(0, column_offset) =
            shape_derivative(master_node, 0) * Mat3::Identity();

        direct_tangent_s.template block<3, 3>(0, column_offset) =
            shape_derivative(master_node, 1) * Mat3::Identity();
    }

    ProjectionMap projection_derivative = ProjectionMap::Zero();
    bool          projection_linearization_valid = false;

    if (projection_info.mode == ProjectionMode::Interior) {
        StaticMatrix<2, 2> projection_hessian;

        projection_hessian(0, 0) =
            tangent_r.dot(tangent_r) - difference.dot(tangent_rr);

        projection_hessian(1, 1) =
            tangent_s.dot(tangent_s) - difference.dot(tangent_ss);

        projection_hessian(0, 1) =
            tangent_r.dot(tangent_s) - difference.dot(tangent_rs);

        projection_hessian(1, 0) =
            projection_hessian(0, 1);

        ProjectionMap projection_rhs;

        projection_rhs.row(0) =
            tangent_r.transpose() * direct_difference +
            difference.transpose() * direct_tangent_r;

        projection_rhs.row(1) =
            tangent_s.transpose() * direct_difference +
            difference.transpose() * direct_tangent_s;

        const Precision scale =
            std::max(projection_hessian.norm(), Precision(1e-16));

        const Precision determinant_tolerance =
            Precision(1e-14) * scale * scale;

        if (projection_hessian.allFinite()) {
            diagnostics.minimum_interior_hessian_ratio =
                std::min(
                    diagnostics.minimum_interior_hessian_ratio,
                    projection_hessian.determinant() / (scale * scale)
                );
        }

        const bool valid_hessian =
            projection_hessian.allFinite() &&
            projection_hessian(0, 0) > Precision(0) &&
            projection_hessian.determinant() > determinant_tolerance;

        if (valid_hessian) {
            const Eigen::LDLT<StaticMatrix<2, 2>> solver(
                projection_hessian
            );

            if (solver.info() == Eigen::Success) {
                const ProjectionMap candidate =
                    solver.solve(projection_rhs);

                if (candidate.allFinite()) {
                    projection_derivative          = candidate;
                    projection_linearization_valid = true;
                }
            }
        }
    } else if (projection_info.mode == ProjectionMode::Edge) {
        const Vec2 edge_direction =
            projection_info.edge_direction;

        const StaticMatrix<N, 1> edge_shape_derivative =
            shape_derivative * edge_direction;

        const Vec3 edge_tangent =
            first * edge_direction;

        const Precision dr = edge_direction(0);
        const Precision ds = edge_direction(1);

        const Vec3 edge_curvature =
            dr * dr * tangent_rr +
            ds * ds * tangent_ss +
            Precision(2) * dr * ds * tangent_rs;

        LocalMap direct_edge_tangent = LocalMap::Zero();

        for (Index master_node = 0; master_node < N; ++master_node) {
            const Index column_offset = 3 * (master_node + 1);

            direct_edge_tangent.template block<3, 3>(0, column_offset) =
                edge_shape_derivative(master_node) * Mat3::Identity();
        }

        const LocalRow edge_rhs =
            edge_tangent.transpose() * direct_difference +
            difference.transpose() * direct_edge_tangent;

        const Precision edge_hessian =
            edge_tangent.squaredNorm() -
            difference.dot(edge_curvature);

        const Precision edge_scale =
            std::max(edge_tangent.squaredNorm(), Precision(1e-16));

        if (std::isfinite(edge_hessian)) {
            diagnostics.minimum_edge_hessian_ratio =
                std::min(
                    diagnostics.minimum_edge_hessian_ratio,
                    edge_hessian / edge_scale
                );
        }

        if (std::isfinite(edge_hessian) &&
            edge_hessian > Precision(1e-12) * edge_scale) {
            const LocalRow eta_derivative =
                edge_rhs / edge_hessian;

            projection_derivative =
                edge_direction * eta_derivative;

            projection_linearization_valid =
                projection_derivative.allFinite();
        }
    }

    /*
     * Very small or nearly edge-aligned master facets can produce a formally
     * finite closest-point derivative with a huge natural-coordinate norm. In
     * that case the contact residual remains well defined, but the Newton
     * tangent overpredicts the motion of the projected point and drives the
     * line search into a non-smooth branch. Falling back to a fixed-local
     * linearization keeps the partner and contact force unchanged while using a
     * more conservative small-sliding tangent for this pair.
     */
    if (projection_linearization_valid &&
        projection_derivative.norm() > maximum_projection_derivative_norm) {
        projection_derivative              = ProjectionMap::Zero();
        projection_linearization_valid     = false;
    }

    switch (projection_info.mode) {
        case ProjectionMode::Interior:
            if (projection_linearization_valid) {
                ++diagnostics.interior_consistent;
            } else {
                ++diagnostics.interior_fallback;
            }
            break;

        case ProjectionMode::Edge:
            if (projection_linearization_valid) {
                ++diagnostics.edge_consistent;
            } else {
                ++diagnostics.edge_fallback;
            }
            break;

        case ProjectionMode::Corner:
            ++diagnostics.corner_direct;
            break;
    }

    diagnostics.maximum_projection_derivative =
        std::max(
            diagnostics.maximum_projection_derivative,
            projection_derivative.norm()
        );

    StaticMatrix<3, 2> surface_tangents;
    surface_tangents.col(0) = tangent_r;
    surface_tangents.col(1) = tangent_s;

    const LocalMap difference_derivative =
        direct_difference -
        surface_tangents * projection_derivative;

    const LocalMap tangent_r_derivative =
        direct_tangent_r +
        tangent_rr * projection_derivative.row(0) +
        tangent_rs * projection_derivative.row(1);

    const LocalMap tangent_s_derivative =
        direct_tangent_s +
        tangent_rs * projection_derivative.row(0) +
        tangent_ss * projection_derivative.row(1);

    const Mat3 normal_projector =
        Mat3::Identity() - normal_base * normal_base.transpose();

    const LocalMap normal_derivative_base =
        normal_projector / surface_jacobian *
        (
            -skew(tangent_s) * tangent_r_derivative +
             skew(tangent_r) * tangent_s_derivative
        );

    const LocalMap normal_derivative =
        flip_normal ? -normal_derivative_base : normal_derivative_base;

    const LocalRow gap_derivative =
        normal.transpose() * difference_derivative +
        difference.transpose() * normal_derivative;

    LocalMatrix contact_tangent = LocalMatrix::Zero();

    for (Index local_node = 0; local_node < local_nodes; ++local_node) {
        Precision weight = Precision(1);
        LocalRow weight_derivative = LocalRow::Zero();

        if (local_node > 0) {
            const Index master_node = local_node - 1;

            weight = -shape(master_node);

            weight_derivative = -(
                shape_derivative(master_node, 0) *
                    projection_derivative.row(0) +
                shape_derivative(master_node, 1) *
                    projection_derivative.row(1)
            );
        }

        const LocalMap force_derivative =
            contact_scale *
            (
                weight * normal *
                    (contact_law.derivative * gap_derivative) +
                contact_law.value * normal * weight_derivative +
                contact_law.value * weight * normal_derivative
            );

        contact_tangent.template block<3, local_dofs>(
            3 * local_node,
            0
        ) = force_derivative;
    }

    logging::error(
        contact_tangent.allFinite(),
        "CONTACT: non-finite local contact tangent"
    );

    const Precision tangent_norm =
        contact_tangent.norm();

    const Precision raw_asymmetry =
        (contact_tangent - contact_tangent.transpose()).norm() /
        std::max(tangent_norm, std::numeric_limits<Precision>::epsilon());

    diagnostics.maximum_raw_asymmetry =
        std::max(
            diagnostics.maximum_raw_asymmetry,
            raw_asymmetry
        );

    diagnostics.maximum_local_tangent_norm =
        std::max(
            diagnostics.maximum_local_tangent_norm,
            tangent_norm
        );

    scatter_contact_tangent<N>(
        local_node_ids,
        contact_tangent,
        system_nodal_dofs,
        triplets
    );
}

void dispatch_contact_pair(
    const model::SurfaceInterface::Ptr& surface,
    ID                                  slave_node_id,
    const Vec2&                         master_local,
    const ProjectionInfo&               projection_info,
    Precision                           penalty,
    Precision                           clearance,
    Precision                           slave_weight,
    bool                                flip_normal,
    SystemDofIds&                       system_nodal_dofs,
    const model::Field&                 node_coords,
    model::NodeData&                    nodal_forces,
    TripletList&                        triplets,
    ContactDiagnostics&                 diagnostics
) {
    switch (surface->n_nodes) {
        case 3: {
            const auto* typed_surface =
                dynamic_cast<const model::Surface<3>*>(surface.get());

            logging::error(
                typed_surface != nullptr,
                "CONTACT: failed to cast three-node master surface"
            );

            assemble_contact_pair<3>(
                *typed_surface,
                slave_node_id,
                master_local,
                projection_info,
                penalty,
                clearance,
                slave_weight,
                flip_normal,
                system_nodal_dofs,
                node_coords,
                nodal_forces,
                triplets,
                diagnostics
            );
            return;
        }

        case 4: {
            const auto* typed_surface =
                dynamic_cast<const model::Surface<4>*>(surface.get());

            logging::error(
                typed_surface != nullptr,
                "CONTACT: failed to cast four-node master surface"
            );

            assemble_contact_pair<4>(
                *typed_surface,
                slave_node_id,
                master_local,
                projection_info,
                penalty,
                clearance,
                slave_weight,
                flip_normal,
                system_nodal_dofs,
                node_coords,
                nodal_forces,
                triplets,
                diagnostics
            );
            return;
        }

        case 6: {
            const auto* typed_surface =
                dynamic_cast<const model::Surface<6>*>(surface.get());

            logging::error(
                typed_surface != nullptr,
                "CONTACT: failed to cast six-node master surface"
            );

            assemble_contact_pair<6>(
                *typed_surface,
                slave_node_id,
                master_local,
                projection_info,
                penalty,
                clearance,
                slave_weight,
                flip_normal,
                system_nodal_dofs,
                node_coords,
                nodal_forces,
                triplets,
                diagnostics
            );
            return;
        }

        case 8: {
            const auto* typed_surface =
                dynamic_cast<const model::Surface<8>*>(surface.get());

            logging::error(
                typed_surface != nullptr,
                "CONTACT: failed to cast eight-node master surface"
            );

            assemble_contact_pair<8>(
                *typed_surface,
                slave_node_id,
                master_local,
                projection_info,
                penalty,
                clearance,
                slave_weight,
                flip_normal,
                system_nodal_dofs,
                node_coords,
                nodal_forces,
                triplets,
                diagnostics
            );
            return;
        }

        default:
            logging::error(
                false,
                "CONTACT: unsupported master surface with ",
                surface->n_nodes,
                " nodes"
            );
    }
}

} // namespace

Contact::Contact(
    model::SurfaceRegion::Ptr master,
    model::NodeRegion::Ptr    slave,
    Precision                 search_distance,
    Precision                 penalty_stiffness,
    Precision                 contact_clearance,
    bool                      flip_master_normal
)
    : master_surfaces(std::move(master)),
      slave_nodes(std::move(slave)),
      slave_surfaces(nullptr),
      distance(search_distance),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(
        distance >= Precision(0),
        "CONTACT: DISTANCE must be non-negative"
    );

    logging::error(
        penalty > Precision(0),
        "CONTACT: PENALTY must be positive"
    );
}

Contact::Contact(
    model::SurfaceRegion::Ptr master,
    model::SurfaceRegion::Ptr slave,
    Precision                 search_distance,
    Precision                 penalty_stiffness,
    Precision                 contact_clearance,
    bool                      flip_master_normal
)
    : master_surfaces(std::move(master)),
      slave_nodes(nullptr),
      slave_surfaces(std::move(slave)),
      distance(search_distance),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(
        distance >= Precision(0),
        "CONTACT: DISTANCE must be non-negative"
    );

    logging::error(
        penalty > Precision(0),
        "CONTACT: PENALTY must be positive"
    );
}

void Contact::begin_trial(
    bool freeze_partners,
    bool freeze_after_update
) const {
    const AssemblyState& parent_state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    runtime_state.trials.push_back({
        parent_state,
        freeze_partners,
        freeze_after_update
    });
}

void Contact::commit_trial() const {
    logging::error(
        !runtime_state.trials.empty(),
        "CONTACT: no active trial state to commit"
    );

    AssemblyState accepted_state =
        std::move(runtime_state.trials.back().state);

    runtime_state.trials.pop_back();

    if (runtime_state.trials.empty()) {
        runtime_state.committed = std::move(accepted_state);
    } else {
        runtime_state.trials.back().state = std::move(accepted_state);
    }
}

void Contact::rollback_trial() const {
    logging::error(
        !runtime_state.trials.empty(),
        "CONTACT: no active trial state to roll back"
    );

    runtime_state.trials.pop_back();
}

bool Contact::partner_signature_changed() const {
    const AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    return state.last_signature_changed;
}

void Contact::assemble(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList&      triplets
) const {
    logging::error(
        model_data.positions != nullptr,
        "CONTACT: positions field not set in model data"
    );

    const model::Field& node_coords =
        *model_data.positions;

    auto& surfaces =
        model_data.surfaces;

    const auto assembly_start =
        std::chrono::steady_clock::now();

    const std::size_t initial_triplet_count =
        triplets.size();

    ContactDiagnostics diagnostics;

    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    // Increments start with one closest-surface search and then freeze the
    // selected master surface. Natural coordinates continue to be reprojected
    // on that surface so the contact point can slide within the element.
    const bool freeze_surface_partners =
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_surface_partners;

    std::unordered_map<ID, ID> current_partners;
    std::unordered_set<ID>     current_active_slaves;

    const Precision search_radius =
        distance + std::max(clearance, Precision(0));

    BvhAabb bvh(search_radius);

    if (!runtime_state.master_topology_initialized) {
        runtime_state.master_edge_neighbors =
            build_master_patch_topology(master_surfaces, surfaces);

        runtime_state.master_topology_initialized = true;
    }

    const auto& master_edge_neighbors =
        runtime_state.master_edge_neighbors;

    for (ID surface_id : *master_surfaces) {
        if (surface_id < 0 ||
            static_cast<std::size_t>(surface_id) >= surfaces.size()) {
            continue;
        }

        const auto& surface =
            surfaces[static_cast<std::size_t>(surface_id)];

        if (!surface) {
            continue;
        }

        bvh.add_element(
            surface_id,
            node_coords,
            surface->nodes(),
            static_cast<int>(surface->n_nodes)
        );
    }

    bvh.finalize();

    if (!bvh.valid()) {
        std::cout
            << "[CONTACT]"
            << " call=" << runtime_state.call
            << " invalid_bvh=1"
            << '\n';
        return;
    }

    if (!runtime_state.slave_nodes_initialized) {
        runtime_state.slave_node_ids =
            collect_slave_nodes(
                slave_nodes,
                slave_surfaces,
                model_data
            );

        runtime_state.slave_nodes_initialized = true;
    }

    const auto& slave_node_ids =
        runtime_state.slave_node_ids;

    diagnostics.slave_nodes =
        static_cast<Index>(slave_node_ids.size());

    if (!runtime_state.slave_weights_initialized) {
        runtime_state.slave_tributary_areas =
            build_slave_tributary_areas(
                slave_surfaces,
                model_data,
                node_coords
            );

        runtime_state.slave_weights_initialized = true;
    }

    const auto& slave_tributary_areas =
        runtime_state.slave_tributary_areas;

    current_partners.reserve(slave_node_ids.size());
    current_active_slaves.reserve(slave_node_ids.size());

    std::vector<ID> candidates;
    candidates.reserve(64);

    for (ID slave_node_id : slave_node_ids) {
        const Vec3 slave_position =
            node_coords.row_vec3(static_cast<Index>(slave_node_id));

        Precision best_distance =
            std::numeric_limits<Precision>::max();

        ID   best_surface_id = ID(-1);
        Vec2 best_local      = Vec2::Zero();

        const auto previous_partner =
            state.partners.find(slave_node_id);

        if (freeze_surface_partners) {
            /*
             * During Newton and line-search assemblies, keep the surface
             * partner selected by the current active-set state. The natural
             * coordinates are still projected on that fixed surface in the
             * current configuration so that the contact point can slide within
             * the element and the tangent contains the local projection
             * derivative.
             */
            if (previous_partner == state.partners.end()) {
                ++diagnostics.no_partner;
                continue;
            }

            best_surface_id = previous_partner->second;

            if (best_surface_id < 0 ||
                static_cast<std::size_t>(best_surface_id) >= surfaces.size() ||
                !surfaces[static_cast<std::size_t>(best_surface_id)]) {
                ++diagnostics.no_partner;
                continue;
            }

            const auto& fixed_surface =
                surfaces[static_cast<std::size_t>(best_surface_id)];

            best_local = fixed_surface->global_to_local(
                slave_position,
                node_coords,
                true
            );

            if (!best_local.allFinite() ||
                !projection_is_inside_closed_domain(
                    fixed_surface->n_nodes,
                    best_local
                )) {
                ++diagnostics.invalid_projections;
                ++diagnostics.no_partner;
                continue;
            }

            ++diagnostics.valid_projections;

            const Vec3 fixed_master_position =
                fixed_surface->local_to_global(best_local, node_coords);

            best_distance =
                (slave_position - fixed_master_position).norm();

            if (!std::isfinite(best_distance)) {
                ++diagnostics.invalid_projections;
                ++diagnostics.no_partner;
                continue;
            }
        } else {
            if (previous_partner != state.partners.end()) {
                const PatchProjection patch_projection =
                    walk_master_patch(
                        previous_partner->second,
                        slave_position,
                        node_coords,
                        surfaces,
                        master_edge_neighbors
                    );

                if (patch_projection.valid &&
                    patch_projection.distance <= search_radius) {
                    best_surface_id = patch_projection.surface_id;
                    best_local      = patch_projection.local;
                    best_distance   = patch_projection.distance;
                }
            }

            if (best_surface_id < 0) {
                const auto& candidate_ids =
                    bvh.query_point(slave_position, &candidates);

                const Index candidate_count =
                    static_cast<Index>(candidate_ids.size());

                diagnostics.candidate_surfaces += candidate_count;
                diagnostics.maximum_candidates = std::max(
                    diagnostics.maximum_candidates,
                    candidate_count
                );

                if (candidate_ids.empty()) {
                    ++diagnostics.zero_candidate_slaves;
                }

                /*
                 * When no previous partner can be tracked through the local
                 * patch, fall back to the bounded closest-point BVH search.
                 * The resulting partner is frozen for the following Newton
                 * and line-search assemblies of the current active-set state.
                 */
                for (ID surface_id : candidate_ids) {
                    if (surface_id < 0 ||
                        static_cast<std::size_t>(surface_id) >= surfaces.size()) {
                        continue;
                    }

                    const auto& surface =
                        surfaces[static_cast<std::size_t>(surface_id)];

                    if (!surface) {
                        continue;
                    }

                    const Vec2 local = surface->global_to_local(
                        slave_position,
                        node_coords,
                        true
                    );

                    if (!local.allFinite() ||
                        !projection_is_inside_closed_domain(
                            surface->n_nodes,
                            local
                        )) {
                        ++diagnostics.invalid_projections;
                        continue;
                    }

                    ++diagnostics.valid_projections;

                    const Vec3 master_position =
                        surface->local_to_global(local, node_coords);

                    const Precision distance_to_surface =
                        (slave_position - master_position).norm();

                    if (!std::isfinite(distance_to_surface) ||
                        distance_to_surface > search_radius) {
                        ++diagnostics.distance_rejections;
                        continue;
                    }

                    constexpr Precision selection_tolerance = Precision(1e-8);

                    const bool first_candidate = best_surface_id < 0;
                    const bool strictly_closer =
                        distance_to_surface < best_distance - selection_tolerance;
                    const bool equal_distance_lower_id =
                        !first_candidate &&
                        std::abs(distance_to_surface - best_distance) <=
                            selection_tolerance &&
                        surface_id < best_surface_id;

                    if (first_candidate ||
                        strictly_closer ||
                        equal_distance_lower_id) {
                        best_distance   = distance_to_surface;
                        best_surface_id = surface_id;
                        best_local      = local;
                    }
                }
            }
        }

        if (best_surface_id < 0) {
            ++diagnostics.no_partner;
            continue;
        }

        /*
         * Faceted master patches often contain two adjacent representations of
         * practically the same closest point. Keep a valid previous partner
         * unless the candidate is clearly closer; this prevents an outer
         * active-set cycle caused by numerically equal edge ownership.
         */
        if (!freeze_surface_partners &&
            previous_partner != state.partners.end() &&
            previous_partner->second != best_surface_id) {
            const Precision switch_tolerance =
                std::max(
                    partner_switch_distance_tolerance,
                    search_radius * Precision(1e-3)
                );

            const PatchProjection previous_projection =
                project_on_surface(
                    previous_partner->second,
                    slave_position,
                    node_coords,
                    surfaces
                );

            if (previous_projection.valid &&
                previous_projection.distance <= search_radius &&
                previous_projection.distance <= best_distance + switch_tolerance) {
                best_surface_id = previous_projection.surface_id;
                best_local      = previous_projection.local;
                best_distance   = previous_projection.distance;
            }
        }

        diagnostics.maximum_closest_distance =
            std::max(
                diagnostics.maximum_closest_distance,
                best_distance
            );

        const auto& best_surface =
            surfaces[static_cast<std::size_t>(best_surface_id)];

        const Vec3 master_position =
            best_surface->local_to_global(best_local, node_coords);

        const ProjectionInfo projection_info =
            classify_projection(
                best_surface->n_nodes,
                best_local
            );

        Vec3 normal =
            best_surface->normal(node_coords, best_local);

        if (flip_normal) {
            normal = -normal;
        }

        logging::error(
            normal.allFinite(),
            "CONTACT: surface ",
            best_surface_id,
            " has invalid normal"
        );

        const Precision gap =
            (slave_position - master_position).dot(normal) -
            clearance;

        const auto slave_area_it =
            slave_tributary_areas.find(slave_node_id);

        const Precision slave_weight =
            slave_area_it != slave_tributary_areas.end()
                ? slave_area_it->second
                : Precision(1);

        const PenaltyLaw contact_law =
            evaluate_penalty_law(gap);

        // Freeze only the selected master surface, not contact activity or
        // natural coordinates. Open pairs remain stored so their gap may become
        // negative in a later Newton or line-search trial.
        current_partners[slave_node_id] = best_surface_id;

        if (gap >= Precision(0)) {
            ++diagnostics.open_closest_partner;
            continue;
        }

        ++diagnostics.active_contacts;
        current_active_slaves.insert(slave_node_id);

        const Precision penetration =
            -gap;

        const Precision slave_force =
            std::abs(
                slave_weight *
                penalty *
                contact_law.value
            );

        diagnostics.slave_force_squared_sum +=
            slave_force * slave_force;

        diagnostics.maximum_slave_force =
            std::max(
                diagnostics.maximum_slave_force,
                slave_force
            );

        if (penetration > diagnostics.maximum_penetration) {
            diagnostics.maximum_penetration = penetration;
            diagnostics.worst_slave         = slave_node_id;
            diagnostics.worst_surface       = best_surface_id;
            diagnostics.worst_local         = best_local;
            diagnostics.worst_mode          = projection_info.mode;
            diagnostics.worst_gap           = gap;
            diagnostics.worst_distance      = best_distance;
        }

        if (state.active_slaves.find(slave_node_id) ==
            state.active_slaves.end()) {
            ++diagnostics.activations;
        }

        if (previous_partner != state.partners.end() &&
            previous_partner->second != best_surface_id) {
            ++diagnostics.partner_switches;

            const ID old_surface_id =
                previous_partner->second;

            if (old_surface_id >= 0 &&
                static_cast<std::size_t>(old_surface_id) < surfaces.size() &&
                surfaces[static_cast<std::size_t>(old_surface_id)]) {
                const auto& old_surface =
                    surfaces[static_cast<std::size_t>(old_surface_id)];

                const Vec2 old_local =
                    old_surface->global_to_local(
                        slave_position,
                        node_coords,
                        true
                    );

                const Vec3 old_master_position =
                    old_surface->local_to_global(
                        old_local,
                        node_coords
                    );

                Vec3 old_normal =
                    old_surface->normal(
                        node_coords,
                        old_local
                    );

                if (flip_normal) {
                    old_normal = -old_normal;
                }

                const Precision old_distance =
                    (slave_position - old_master_position).norm();

                const Precision old_gap =
                    (slave_position - old_master_position).dot(old_normal) -
                    clearance;

                const ProjectionInfo old_projection_info =
                    classify_projection(
                        old_surface->n_nodes,
                        old_local
                    );

                const auto previous_flags =
                    std::cout.flags();

                const auto previous_precision =
                    std::cout.precision();

                std::cout
                    << std::scientific
                    << std::setprecision(9)
                    << "[CONTACT_SWITCH]"
                    << " call="       << runtime_state.call
                    << " slave="      << slave_node_id
                    << " old="        << old_surface_id
                    << " new="        << best_surface_id
                    << " old_mode="   << projection_mode_name(old_projection_info.mode)
                    << " new_mode="   << projection_mode_name(projection_info.mode)
                    << " old_dist="   << old_distance
                    << " new_dist="   << best_distance
                    << " delta_dist=" << best_distance - old_distance
                    << " old_gap="    << old_gap
                    << " new_gap="    << gap
                    << " old_local=(" << old_local(0)
                    << ','             << old_local(1)
                    << ")"
                    << " new_local=(" << best_local(0)
                    << ','             << best_local(1)
                    << ")"
                    << '\n';

                std::cout.flags(previous_flags);
                std::cout.precision(previous_precision);
            }
        }

        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(slave_node_id)
            )
        );

        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(best_surface_id)
            )
        );

        dispatch_contact_pair(
            best_surface,
            slave_node_id,
            best_local,
            projection_info,
            penalty,
            clearance,
            slave_weight,
            flip_normal,
            system_nodal_dofs,
            node_coords,
            nodal_forces,
            triplets,
            diagnostics
        );
    }

    for (ID slave_node_id : state.active_slaves) {
        if (current_active_slaves.find(slave_node_id) ==
            current_active_slaves.end()) {
            ++diagnostics.deactivations;
        }
    }

    const bool signature_changed =
        state.previous_signature != 0 &&
        state.previous_signature != diagnostics.signature;

    const auto active_change =
        static_cast<std::int64_t>(diagnostics.active_contacts) -
        static_cast<std::int64_t>(state.previous_active);

    state.partners           = std::move(current_partners);
    state.active_slaves      = std::move(current_active_slaves);
    state.previous_signature = diagnostics.signature;
    state.previous_active    = diagnostics.active_contacts;
    state.last_signature_changed = signature_changed;

    if (!freeze_surface_partners &&
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_after_update) {
        runtime_state.trials.back().freeze_surface_partners = true;
        runtime_state.trials.back().freeze_after_update     = false;
    }

    const Precision contact_force_norm =
        std::sqrt(diagnostics.slave_force_squared_sum);

    const auto assembly_end =
        std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            assembly_end - assembly_start
        ).count();

    const std::size_t added_triplets =
        triplets.size() - initial_triplet_count;

    const Precision average_candidates =
        diagnostics.slave_nodes > 0
            ? static_cast<Precision>(diagnostics.candidate_surfaces) /
              static_cast<Precision>(diagnostics.slave_nodes)
            : Precision(0);

    const auto previous_flags =
        std::cout.flags();

    const auto previous_precision =
        std::cout.precision();

    std::cout
        << std::scientific
        << std::setprecision(3)
        << "[CONTACT]"
        << " call="          << runtime_state.call
        << " frozen="        << (freeze_surface_partners ? 1 : 0)
        << " active="        << diagnostics.active_contacts
        << " d_active="      << active_change
        << " activated="     << diagnostics.activations
        << " deactivated="   << diagnostics.deactivations
        << " switches="      << diagnostics.partner_switches
        << " changed="       << (signature_changed ? 1 : 0)
        << " interior="      << diagnostics.interior_consistent
        << " interior_fb="   << diagnostics.interior_fallback
        << " edge="          << diagnostics.edge_consistent
        << " edge_fb="       << diagnostics.edge_fallback
        << " corner="        << diagnostics.corner_direct
        << " no_partner="    << diagnostics.no_partner
        << " zero_cand="     << diagnostics.zero_candidate_slaves
        << " invalid_proj="  << diagnostics.invalid_projections
        << " distance_rej="  << diagnostics.distance_rejections
        << " open="          << diagnostics.open_closest_partner
        << " cand_avg="      << average_candidates
        << " cand_max="      << diagnostics.maximum_candidates
        << " max_dist="      << diagnostics.maximum_closest_distance
        << " max_pen="       << diagnostics.maximum_penetration
        << " max_force="     << diagnostics.maximum_slave_force
        << " force_norm="    << contact_force_norm
        << " min_Js="        << diagnostics.minimum_surface_jacobian
        << " min_H2="        << diagnostics.minimum_interior_hessian_ratio
        << " min_H1="        << diagnostics.minimum_edge_hessian_ratio
        << " max_dproj="     << diagnostics.maximum_projection_derivative
        << " max_asym="      << diagnostics.maximum_raw_asymmetry
        << " max_Kc="        << diagnostics.maximum_local_tangent_norm
        << " triplets="      << added_triplets
        << " ms="            << elapsed_ms
        << " signature="     << diagnostics.signature
        << '\n';

    if (diagnostics.active_contacts > 0) {
        std::cout
            << std::scientific
            << std::setprecision(6)
            << "[CONTACT_MAX]"
            << " call="     << runtime_state.call
            << " slave="    << diagnostics.worst_slave
            << " surface="  << diagnostics.worst_surface
            << " mode="     << projection_mode_name(diagnostics.worst_mode)
            << " local=("   << diagnostics.worst_local(0)
            << ","          << diagnostics.worst_local(1)
            << ")"
            << " gap="      << diagnostics.worst_gap
            << " pen="      << diagnostics.maximum_penetration
            << " distance=" << diagnostics.worst_distance
            << '\n';
    }

    std::cout.flags(previous_flags);
    std::cout.precision(previous_precision);
}

} // namespace constraint