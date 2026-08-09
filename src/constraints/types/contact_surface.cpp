/**
 * @file contact_surface.cpp
 * @brief Implements surface-integration-point slave contact.
 *
 * SurfaceRegion slaves are integrated with the native quadrature rule of each
 * slave surface. Every quadrature point acts as one transactional contact point:
 * its current slave position is projected onto the master patch, the normal
 * augmented-Lagrange law is evaluated there, and the resulting residual and
 * tangent are distributed consistently through both slave and master shape
 * functions.
 *
 * The integration measure is evaluated in the reference slave geometry. This
 * preserves the fixed tributary weighting used by the nodal formulation while
 * replacing nodal lumping by an actual surface quadrature. Explicit NodeRegion
 * slaves continue to use Contact::assemble() in contact.cpp.
 *
 * @author Finn Eggers
 * @date 09.08.2026
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
constexpr Precision projection_selection_tolerance     = Precision(1e-8);
constexpr bool      print_surface_contact_summary       = true;

enum class ProjectionMode {
    Interior,
    Edge,
    Corner
};

struct ProjectionInfo {
    ProjectionMode mode           = ProjectionMode::Interior;
    Vec2           edge_direction = Vec2::Zero();
};

struct PatchProjection {
    bool      valid      = false;
    ID        surface_id = ID(-1);
    Vec2      local      = Vec2::Zero();
    Precision distance   = std::numeric_limits<Precision>::max();
    Precision gap        = std::numeric_limits<Precision>::max();
};

struct SurfaceContactDiagnostics {
    Index slave_points             = 0;
    Index candidate_surfaces       = 0;
    Index maximum_candidates       = 0;
    Index zero_candidate_points    = 0;
    Index valid_projections        = 0;
    Index invalid_projections      = 0;
    Index distance_rejections      = 0;
    Index no_partner               = 0;
    Index open_closest_partner     = 0;
    Index active_contacts          = 0;
    Index activations              = 0;
    Index deactivations            = 0;
    Index partner_switches         = 0;
    Index self_contact_rejections  = 0;
    Index active_partner_losses    = 0;
    Index consistent_tangents      = 0;
    Index stabilized_tangents      = 0;

    Precision maximum_penetration       = Precision(0);
    Precision maximum_closest_distance  = Precision(0);
    Precision maximum_point_force       = Precision(0);
    Precision point_force_squared_sum   = Precision(0);
    Precision minimum_master_jacobian   = std::numeric_limits<Precision>::infinity();
    Precision maximum_raw_asymmetry     = Precision(0);
    Precision maximum_local_tangent_norm = Precision(0);

    ID        worst_point   = ID(-1);
    ID        worst_surface = ID(-1);
    Vec2      worst_local   = Vec2::Zero();
    Precision worst_gap     = Precision(0);
    Precision worst_distance = Precision(0);

    std::uint64_t signature = 1469598103934665603ull;
};

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

void hash_value(std::uint64_t& signature, std::uint64_t value) {
    constexpr std::uint64_t prime = 1099511628211ull;

    for (Index byte = 0; byte < 8; ++byte) {
        signature ^= value & 0xffull;
        signature *= prime;
        value >>= 8;
    }
}

bool projection_is_inside_closed_domain(Index number_of_nodes, const Vec2& local) {
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

ProjectionInfo classify_projection(Index number_of_nodes, const Vec2& local) {
    constexpr Precision edge_tolerance = Precision(1e-7);

    ProjectionInfo info;

    if (number_of_nodes == 3 || number_of_nodes == 6) {
        const Precision r = local(0);
        const Precision s = local(1);
        const Precision t = Precision(1) - r - s;

        const bool on_r_edge = std::abs(r) <= edge_tolerance;
        const bool on_s_edge = std::abs(s) <= edge_tolerance;
        const bool on_t_edge = std::abs(t) <= edge_tolerance;

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

        return info;
    }

    if (number_of_nodes == 4 || number_of_nodes == 8) {
        const Precision r = local(0);
        const Precision s = local(1);

        const bool on_r_min_edge = std::abs(r + Precision(1)) <= edge_tolerance;
        const bool on_r_max_edge = std::abs(r - Precision(1)) <= edge_tolerance;
        const bool on_s_min_edge = std::abs(s + Precision(1)) <= edge_tolerance;
        const bool on_s_max_edge = std::abs(s - Precision(1)) <= edge_tolerance;

        const bool on_r_edge = on_r_min_edge || on_r_max_edge;
        const bool on_s_edge = on_s_min_edge || on_s_max_edge;

        if (on_r_edge && on_s_edge) {
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

        return info;
    }

    info.mode = ProjectionMode::Corner;
    return info;
}

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

EdgeKey make_edge_key(const model::SurfaceInterface& surface, Index edge) {
    const auto local_nodes = edge_local_nodes(surface, edge);

    ID first  = surface.nodes()[local_nodes[0]];
    ID second = surface.nodes()[local_nodes[1]];

    const ID middle =
        local_nodes[2] >= 0
            ? surface.nodes()[local_nodes[2]]
            : ID(-1);

    if (second < first) {
        std::swap(first, second);
    }

    return {first, second, middle};
}

ID crossed_edge_index(const model::SurfaceInterface& surface, const Vec2& local) {
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

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
        if (!surface) {
            continue;
        }

        for (Index edge = 0; edge < surface->n_edges; ++edge) {
            const EdgeKey key = make_edge_key(*surface, edge);
            const auto owner = edge_owner.find(key);

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
                          [static_cast<std::size_t>(edge)] = neighbor_id;

            edge_neighbors[static_cast<std::size_t>(neighbor_id)]
                          [static_cast<std::size_t>(neighbor_edge)] = surface_id;
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

bool surfaces_share_node(
    const model::SurfaceInterface& first,
    const model::SurfaceInterface& second
) {
    for (Index i = 0; i < first.n_nodes; ++i) {
        for (Index j = 0; j < second.n_nodes; ++j) {
            if (first.nodes()[i] == second.nodes()[j]) {
                return true;
            }
        }
    }

    return false;
}

PatchProjection evaluate_projection(
    ID                                                surface_id,
    const Vec2&                                       local,
    const Vec3&                                       slave_position,
    const model::Field&                               node_coords,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    Precision                                         clearance,
    bool                                              flip_normal
) {
    PatchProjection projection;

    if (!valid_surface_id(surface_id, surfaces)) {
        return projection;
    }

    const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

    if (!local.allFinite() ||
        !projection_is_inside_closed_domain(surface->n_nodes, local)) {
        return projection;
    }

    const Vec3 master_position = surface->local_to_global(local, node_coords);
    Vec3       normal          = surface->normal(node_coords, local);

    if (flip_normal) {
        normal = -normal;
    }

    const Vec3 difference = slave_position - master_position;

    if (!master_position.allFinite() ||
        !normal.allFinite() ||
        !difference.allFinite()) {
        return projection;
    }

    const Precision distance_to_surface = difference.norm();
    const Precision gap                 = difference.dot(normal) - clearance;

    if (!std::isfinite(distance_to_surface) || !std::isfinite(gap)) {
        return projection;
    }

    projection.valid      = true;
    projection.surface_id = surface_id;
    projection.local      = local;
    projection.distance   = distance_to_surface;
    projection.gap        = gap;

    return projection;
}

PatchProjection project_on_surface(
    ID                                                surface_id,
    const Vec3&                                       slave_position,
    const model::Field&                               node_coords,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    Precision                                         clearance,
    bool                                              flip_normal
) {
    if (!valid_surface_id(surface_id, surfaces)) {
        return {};
    }

    const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
    const Vec2 local = surface->global_to_local(slave_position, node_coords, true);

    return evaluate_projection(
        surface_id,
        local,
        slave_position,
        node_coords,
        surfaces,
        clearance,
        flip_normal
    );
}

bool projection_is_better(
    const PatchProjection& candidate,
    const PatchProjection& current
) {
    if (!candidate.valid) {
        return false;
    }

    if (!current.valid) {
        return true;
    }

    if (candidate.distance < current.distance - projection_selection_tolerance) {
        return true;
    }

    return std::abs(candidate.distance - current.distance) <=
               projection_selection_tolerance &&
           candidate.surface_id < current.surface_id;
}

PatchProjection walk_master_patch(
    ID                                                start_surface_id,
    const Vec3&                                       slave_position,
    const model::Field&                               node_coords,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const std::vector<std::array<ID, 4>>&             edge_neighbors,
    Precision                                         clearance,
    bool                                              flip_normal
) {
    ID current_surface_id = start_surface_id;

    std::unordered_set<ID> visited;
    visited.reserve(
        std::min<std::size_t>(edge_neighbors.size(), std::size_t(32))
    );

    const Index maximum_steps =
        static_cast<Index>(
            std::min<std::size_t>(edge_neighbors.size(), std::size_t(64))
        );

    for (Index step = 0; step < maximum_steps; ++step) {
        if (!valid_surface_id(current_surface_id, surfaces) ||
            visited.find(current_surface_id) != visited.end()) {
            break;
        }

        visited.insert(current_surface_id);

        const auto& surface = surfaces[static_cast<std::size_t>(current_surface_id)];
        const Vec2 unclipped_local =
            surface->global_to_local(slave_position, node_coords, false);

        if (unclipped_local.allFinite() &&
            projection_is_inside_closed_domain(surface->n_nodes, unclipped_local)) {
            return evaluate_projection(
                current_surface_id,
                unclipped_local,
                slave_position,
                node_coords,
                surfaces,
                clearance,
                flip_normal
            );
        }

        const ID crossed_edge =
            unclipped_local.allFinite()
                ? crossed_edge_index(*surface, unclipped_local)
                : ID(-1);

        if (crossed_edge < 0 ||
            static_cast<std::size_t>(current_surface_id) >= edge_neighbors.size()) {
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

    return {};
}

struct ContactLaw {
    Precision value      = Precision(0);
    Precision derivative = Precision(0);
    Precision pressure   = Precision(0);
    bool      active     = false;
};

ContactLaw evaluate_augmented_lagrange_law(
    Precision gap,
    Precision normal_multiplier,
    Precision penalty
) {
    const Precision shifted_gap = gap - normal_multiplier / penalty;

    if (shifted_gap >= Precision(0)) {
        return {};
    }

    return {
        shifted_gap,
        Precision(1),
        -penalty * shifted_gap,
        true
    };
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

template<Index NS, Index NM>
void scatter_surface_contact_tangent(
    const std::array<ID, NS + NM>&                        node_ids,
    const StaticMatrix<3 * (NS + NM), 3 * (NS + NM)>&   tangent,
    SystemDofIds&                                        system_nodal_dofs,
    TripletList&                                         triplets
) {
    constexpr Index local_nodes = NS + NM;
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
        const int global_row = global_dofs[static_cast<std::size_t>(local_row)];
        if (global_row < 0) {
            continue;
        }

        for (Index local_col = 0; local_col < local_dofs; ++local_col) {
            const int global_col = global_dofs[static_cast<std::size_t>(local_col)];
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

template<Index NS, Index NM>
void assemble_surface_contact_pair(
    const model::Surface<NS>& slave_surface,
    const model::Surface<NM>& master_surface,
    const Vec2&               slave_local,
    const Vec2&               master_local,
    const ProjectionInfo&     projection_info,
    Precision                 penalty,
    Precision                 clearance,
    Precision                 normal_multiplier,
    Precision                 integration_weight,
    bool                      flip_normal,
    SystemDofIds&             system_nodal_dofs,
    const model::Field&       node_coords,
    model::NodeData&          nodal_forces,
    TripletList&              triplets,
    SurfaceContactDiagnostics& diagnostics
) {
    constexpr Index local_nodes = NS + NM;
    constexpr Index local_dofs  = 3 * local_nodes;

    using LocalMatrix   = StaticMatrix<local_dofs, local_dofs>;
    using LocalRow      = StaticMatrix<1, local_dofs>;
    using LocalMap      = StaticMatrix<3, local_dofs>;
    using ProjectionMap = StaticMatrix<2, local_dofs>;

    const StaticMatrix<NS, 1> slave_shape =
        slave_surface.shape_function(slave_local(0), slave_local(1));

    const StaticMatrix<NM, 1> master_shape =
        master_surface.shape_function(master_local(0), master_local(1));

    const StaticMatrix<NM, 2> master_shape_derivative =
        master_surface.shape_derivative(master_local(0), master_local(1));

    const StaticMatrix<NM, 3> master_shape_second_derivative =
        master_surface.shape_second_derivative(master_local(0), master_local(1));

    const StaticMatrix<NS, 3> slave_coordinates =
        slave_surface.node_coords_global(node_coords);

    const StaticMatrix<NM, 3> master_coordinates =
        master_surface.node_coords_global(node_coords);

    const Vec3 slave_position =
        slave_coordinates.transpose() * slave_shape;

    const Vec3 master_position =
        master_coordinates.transpose() * master_shape;

    const Vec3 difference = slave_position - master_position;

    const StaticMatrix<3, 2> first =
        master_coordinates.transpose() * master_shape_derivative;

    const Vec3 tangent_r = first.col(0);
    const Vec3 tangent_s = first.col(1);

    const StaticMatrix<3, 3> second =
        master_coordinates.transpose() * master_shape_second_derivative;

    const Vec3 tangent_rr = second.col(0);
    const Vec3 tangent_ss = second.col(1);
    const Vec3 tangent_rs = second.col(2);

    const Vec3 normal_unnormalized = tangent_r.cross(tangent_s);
    const Precision master_jacobian = normal_unnormalized.norm();

    diagnostics.minimum_master_jacobian =
        std::min(diagnostics.minimum_master_jacobian, master_jacobian);

    logging::error(
        std::isfinite(master_jacobian) && master_jacobian > Precision(1e-14),
        "CONTACT: singular master-surface Jacobian"
    );

    const Vec3 normal_base = normal_unnormalized / master_jacobian;
    const Vec3 normal      = flip_normal ? -normal_base : normal_base;
    const Precision gap    = difference.dot(normal) - clearance;

    const ContactLaw contact_law =
        evaluate_augmented_lagrange_law(gap, normal_multiplier, penalty);

    if (!contact_law.active) {
        return;
    }

    std::array<ID, local_nodes> local_node_ids{};

    for (Index slave_node = 0; slave_node < NS; ++slave_node) {
        local_node_ids[static_cast<std::size_t>(slave_node)] =
            slave_surface.nodes()[slave_node];
    }

    for (Index master_node = 0; master_node < NM; ++master_node) {
        local_node_ids[static_cast<std::size_t>(NS + master_node)] =
            master_surface.nodes()[master_node];
    }

    const Precision contact_scale = integration_weight * penalty;
    const Vec3 contact_force = contact_scale * contact_law.value * normal;

    for (Index slave_node = 0; slave_node < NS; ++slave_node) {
        add_translational_force(
            nodal_forces,
            slave_surface.nodes()[slave_node],
            slave_shape(slave_node) * contact_force
        );
    }

    for (Index master_node = 0; master_node < NM; ++master_node) {
        add_translational_force(
            nodal_forces,
            master_surface.nodes()[master_node],
            -master_shape(master_node) * contact_force
        );
    }

    LocalMap direct_difference = LocalMap::Zero();

    for (Index slave_node = 0; slave_node < NS; ++slave_node) {
        direct_difference.template block<3, 3>(0, 3 * slave_node) =
            slave_shape(slave_node) * Mat3::Identity();
    }

    LocalMap direct_tangent_r = LocalMap::Zero();
    LocalMap direct_tangent_s = LocalMap::Zero();

    for (Index master_node = 0; master_node < NM; ++master_node) {
        const Index column_offset = 3 * (NS + master_node);

        direct_difference.template block<3, 3>(0, column_offset) =
            -master_shape(master_node) * Mat3::Identity();

        direct_tangent_r.template block<3, 3>(0, column_offset) =
            master_shape_derivative(master_node, 0) * Mat3::Identity();

        direct_tangent_s.template block<3, 3>(0, column_offset) =
            master_shape_derivative(master_node, 1) * Mat3::Identity();
    }

    ProjectionMap projection_derivative = ProjectionMap::Zero();
    bool projection_linearization_valid = false;

    if (projection_info.mode == ProjectionMode::Interior) {
        StaticMatrix<2, 2> projection_hessian;

        projection_hessian(0, 0) =
            tangent_r.dot(tangent_r) - difference.dot(tangent_rr);

        projection_hessian(1, 1) =
            tangent_s.dot(tangent_s) - difference.dot(tangent_ss);

        projection_hessian(0, 1) =
            tangent_r.dot(tangent_s) - difference.dot(tangent_rs);

        projection_hessian(1, 0) = projection_hessian(0, 1);

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

        const bool valid_hessian =
            projection_hessian.allFinite() &&
            projection_hessian(0, 0) > Precision(0) &&
            projection_hessian.determinant() > determinant_tolerance;

        if (valid_hessian) {
            const Eigen::LDLT<StaticMatrix<2, 2>> solver(projection_hessian);

            if (solver.info() == Eigen::Success) {
                const ProjectionMap candidate = solver.solve(projection_rhs);

                if (candidate.allFinite()) {
                    projection_derivative          = candidate;
                    projection_linearization_valid = true;
                }
            }
        }
    } else if (projection_info.mode == ProjectionMode::Edge) {
        const Vec2 edge_direction = projection_info.edge_direction;

        const StaticMatrix<NM, 1> edge_shape_derivative =
            master_shape_derivative * edge_direction;

        const Vec3 edge_tangent = first * edge_direction;

        const Precision dr = edge_direction(0);
        const Precision ds = edge_direction(1);

        const Vec3 edge_curvature =
            dr * dr * tangent_rr +
            ds * ds * tangent_ss +
            Precision(2) * dr * ds * tangent_rs;

        LocalMap direct_edge_tangent = LocalMap::Zero();

        for (Index master_node = 0; master_node < NM; ++master_node) {
            const Index column_offset = 3 * (NS + master_node);

            direct_edge_tangent.template block<3, 3>(0, column_offset) =
                edge_shape_derivative(master_node) * Mat3::Identity();
        }

        const LocalRow edge_rhs =
            edge_tangent.transpose() * direct_difference +
            difference.transpose() * direct_edge_tangent;

        const Precision edge_hessian =
            edge_tangent.squaredNorm() - difference.dot(edge_curvature);

        const Precision edge_scale =
            std::max(edge_tangent.squaredNorm(), Precision(1e-16));

        if (std::isfinite(edge_hessian) &&
            edge_hessian > Precision(1e-12) * edge_scale) {
            const LocalRow eta_derivative = edge_rhs / edge_hessian;
            projection_derivative = edge_direction * eta_derivative;
            projection_linearization_valid = projection_derivative.allFinite();
        }
    }

    if (projection_linearization_valid &&
        projection_derivative.norm() > maximum_projection_derivative_norm) {
        projection_derivative          = ProjectionMap::Zero();
        projection_linearization_valid = false;
    }

    LocalMatrix contact_tangent = LocalMatrix::Zero();

    const bool use_consistent_tangent =
        projection_info.mode == ProjectionMode::Interior &&
        projection_linearization_valid;

    if (use_consistent_tangent) {
        StaticMatrix<3, 2> surface_tangents;
        surface_tangents.col(0) = tangent_r;
        surface_tangents.col(1) = tangent_s;

        const LocalMap difference_derivative =
            direct_difference - surface_tangents * projection_derivative;

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
            normal_projector / master_jacobian *
            (
                -skew(tangent_s) * tangent_r_derivative +
                 skew(tangent_r) * tangent_s_derivative
            );

        const LocalMap normal_derivative =
            flip_normal ? -normal_derivative_base : normal_derivative_base;

        const LocalRow gap_derivative =
            normal.transpose() * difference_derivative +
            difference.transpose() * normal_derivative;

        for (Index local_node = 0; local_node < local_nodes; ++local_node) {
            Precision weight = Precision(0);
            LocalRow weight_derivative = LocalRow::Zero();

            if (local_node < NS) {
                weight = slave_shape(local_node);
            } else {
                const Index master_node = local_node - NS;
                weight = -master_shape(master_node);

                weight_derivative = -(
                    master_shape_derivative(master_node, 0) *
                        projection_derivative.row(0) +
                    master_shape_derivative(master_node, 1) *
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

        ++diagnostics.consistent_tangents;
    } else {
        std::array<Precision, local_nodes> weights{};

        for (Index slave_node = 0; slave_node < NS; ++slave_node) {
            weights[static_cast<std::size_t>(slave_node)] =
                slave_shape(slave_node);
        }

        for (Index master_node = 0; master_node < NM; ++master_node) {
            weights[static_cast<std::size_t>(NS + master_node)] =
                -master_shape(master_node);
        }

        const Mat3 normal_stiffness =
            contact_scale * contact_law.derivative *
            normal * normal.transpose();

        for (Index local_row = 0; local_row < local_nodes; ++local_row) {
            for (Index local_col = 0; local_col < local_nodes; ++local_col) {
                contact_tangent.template block<3, 3>(
                    3 * local_row,
                    3 * local_col
                ) = weights[static_cast<std::size_t>(local_row)] *
                    weights[static_cast<std::size_t>(local_col)] *
                    normal_stiffness;
            }
        }

        ++diagnostics.stabilized_tangents;
    }

    logging::error(
        contact_tangent.allFinite(),
        "CONTACT: non-finite surface-integration-point tangent"
    );

    const Precision tangent_norm = contact_tangent.norm();
    const Precision raw_asymmetry =
        (contact_tangent - contact_tangent.transpose()).norm() /
        std::max(tangent_norm, std::numeric_limits<Precision>::epsilon());

    diagnostics.maximum_raw_asymmetry =
        std::max(diagnostics.maximum_raw_asymmetry, raw_asymmetry);

    diagnostics.maximum_local_tangent_norm =
        std::max(diagnostics.maximum_local_tangent_norm, tangent_norm);

    scatter_surface_contact_tangent<NS, NM>(
        local_node_ids,
        contact_tangent,
        system_nodal_dofs,
        triplets
    );
}

template<Index NS>
void dispatch_surface_contact_pair(
    const model::Surface<NS>&             slave_surface,
    const model::SurfaceInterface::Ptr&   master_surface,
    const Vec2&                           slave_local,
    const Vec2&                           master_local,
    const ProjectionInfo&                 projection_info,
    Precision                             penalty,
    Precision                             clearance,
    Precision                             normal_multiplier,
    Precision                             integration_weight,
    bool                                  flip_normal,
    SystemDofIds&                         system_nodal_dofs,
    const model::Field&                   node_coords,
    model::NodeData&                      nodal_forces,
    TripletList&                          triplets,
    SurfaceContactDiagnostics&            diagnostics
) {
    switch (master_surface->n_nodes) {
        case 3: {
            const auto* typed =
                dynamic_cast<const model::Surface<3>*>(master_surface.get());
            logging::error(typed != nullptr,
                "CONTACT: failed to cast three-node master surface");
            assemble_surface_contact_pair<NS, 3>(
                slave_surface, *typed, slave_local, master_local,
                projection_info, penalty, clearance, normal_multiplier,
                integration_weight, flip_normal, system_nodal_dofs,
                node_coords, nodal_forces, triplets, diagnostics
            );
            return;
        }

        case 4: {
            const auto* typed =
                dynamic_cast<const model::Surface<4>*>(master_surface.get());
            logging::error(typed != nullptr,
                "CONTACT: failed to cast four-node master surface");
            assemble_surface_contact_pair<NS, 4>(
                slave_surface, *typed, slave_local, master_local,
                projection_info, penalty, clearance, normal_multiplier,
                integration_weight, flip_normal, system_nodal_dofs,
                node_coords, nodal_forces, triplets, diagnostics
            );
            return;
        }

        case 6: {
            const auto* typed =
                dynamic_cast<const model::Surface<6>*>(master_surface.get());
            logging::error(typed != nullptr,
                "CONTACT: failed to cast six-node master surface");
            assemble_surface_contact_pair<NS, 6>(
                slave_surface, *typed, slave_local, master_local,
                projection_info, penalty, clearance, normal_multiplier,
                integration_weight, flip_normal, system_nodal_dofs,
                node_coords, nodal_forces, triplets, diagnostics
            );
            return;
        }

        case 8: {
            const auto* typed =
                dynamic_cast<const model::Surface<8>*>(master_surface.get());
            logging::error(typed != nullptr,
                "CONTACT: failed to cast eight-node master surface");
            assemble_surface_contact_pair<NS, 8>(
                slave_surface, *typed, slave_local, master_local,
                projection_info, penalty, clearance, normal_multiplier,
                integration_weight, flip_normal, system_nodal_dofs,
                node_coords, nodal_forces, triplets, diagnostics
            );
            return;
        }

        default:
            logging::error(
                false,
                "CONTACT: unsupported master surface with ",
                master_surface->n_nodes,
                " nodes"
            );
    }
}

} // namespace

void Contact::assemble_surface(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList&      triplets
) const {
    logging::error(
        slave_surfaces != nullptr,
        "CONTACT: surface assembly requires a SurfaceRegion slave"
    );

    logging::error(
        slave_nodes == nullptr,
        "CONTACT: surface assembly cannot be used with a NodeRegion slave"
    );

    logging::error(
        model_data.positions != nullptr,
        "CONTACT: positions field not set in model data"
    );

    logging::error(
        model_data.positions_reference != nullptr,
        "CONTACT: reference positions field not set for surface integration"
    );

    const model::Field& node_coords = *model_data.positions;
    const model::Field& reference_coords = *model_data.positions_reference;
    auto& surfaces = model_data.surfaces;

    const auto assembly_start = std::chrono::steady_clock::now();
    const std::size_t initial_triplet_count = triplets.size();

    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    const bool freeze_surface_partners =
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_surface_partners;

    const Precision search_radius =
        distance + std::max(clearance, Precision(0));

    BvhAabb bvh(search_radius);

    if (!runtime_state.master_topology_initialized) {
        runtime_state.master_edge_neighbors =
            build_master_patch_topology(master_surfaces, surfaces);
        runtime_state.master_topology_initialized = true;
    }

    const auto& master_edge_neighbors = runtime_state.master_edge_neighbors;

    for (ID surface_id : *master_surfaces) {
        if (!valid_surface_id(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

        bvh.add_element(
            surface_id,
            node_coords,
            surface->nodes(),
            static_cast<int>(surface->n_nodes)
        );
    }

    bvh.finalize();

    std::unordered_map<ID, ID> current_partners;
    std::unordered_set<ID>     current_active_slaves;
    std::unordered_map<ID, Precision> current_multipliers;
    std::unordered_map<ID, Precision> current_gaps;
    std::unordered_map<ID, Precision> current_characteristic_lengths;

    std::vector<ID> candidates;
    candidates.reserve(64);

    SurfaceContactDiagnostics diagnostics;

    std::vector<ID> slave_surface_ids;
    slave_surface_ids.reserve(static_cast<std::size_t>(slave_surfaces->size()));

    for (ID surface_id : *slave_surfaces) {
        slave_surface_ids.push_back(surface_id);
    }

    std::sort(slave_surface_ids.begin(), slave_surface_ids.end());
    slave_surface_ids.erase(
        std::unique(slave_surface_ids.begin(), slave_surface_ids.end()),
        slave_surface_ids.end()
    );

    ID next_point_id = 0;

    auto process_slave_surface = [&](const auto& slave_surface, ID slave_surface_id) {
        using SlaveSurface = std::decay_t<decltype(slave_surface)>;
        constexpr Index NS = SlaveSurface::num_nodes;

        const auto current_slave_coordinates =
            slave_surface.node_coords_global(node_coords);

        const auto reference_slave_coordinates =
            slave_surface.node_coords_global(reference_coords);

        const auto& scheme = slave_surface.integration_scheme();

        for (Index local_ip = 0; local_ip < scheme.count(); ++local_ip) {
            const auto point = scheme.get_point(local_ip);
            const ID point_id = next_point_id++;

            ++diagnostics.slave_points;

            const Vec2 slave_local(point.r, point.s);
            const StaticMatrix<NS, 1> slave_shape =
                slave_surface.shape_function(point.r, point.s);

            const Vec3 slave_position =
                current_slave_coordinates.transpose() * slave_shape;

            const auto slave_jacobian =
                slave_surface.jacobian(
                    reference_slave_coordinates,
                    point.r,
                    point.s
                );

            const Precision integration_weight =
                slave_jacobian.col(0).cross(slave_jacobian.col(1)).norm() *
                point.w;

            logging::error(
                std::isfinite(integration_weight) &&
                integration_weight > Precision(0),
                "CONTACT: invalid slave integration weight on surface ",
                slave_surface_id,
                ", integration point ",
                local_ip
            );

            const auto previous_partner = state.partners.find(point_id);
            const auto previous_multiplier =
                state.normal_multipliers.find(point_id);

            const Precision normal_multiplier =
                previous_multiplier != state.normal_multipliers.end()
                    ? std::max(previous_multiplier->second, Precision(0))
                    : Precision(0);

            const bool previously_active =
                state.active_slaves.find(point_id) != state.active_slaves.end();

            const bool track_established_partner =
                previous_partner != state.partners.end() &&
                (previously_active || normal_multiplier > Precision(0));

            PatchProjection best_projection;

            if (freeze_surface_partners) {
                if (previous_partner != state.partners.end() &&
                    valid_surface_id(previous_partner->second, surfaces)) {
                    const auto& fixed_surface =
                        surfaces[static_cast<std::size_t>(previous_partner->second)];

                    if (surfaces_share_node(slave_surface, *fixed_surface)) {
                        ++diagnostics.self_contact_rejections;
                    } else {
                        best_projection = project_on_surface(
                            previous_partner->second,
                            slave_position,
                            node_coords,
                            surfaces,
                            clearance,
                            flip_normal
                        );

                        if (best_projection.valid) {
                            ++diagnostics.valid_projections;
                        } else {
                            ++diagnostics.invalid_projections;
                        }
                    }
                }
            } else if (track_established_partner) {
                best_projection = walk_master_patch(
                    previous_partner->second,
                    slave_position,
                    node_coords,
                    surfaces,
                    master_edge_neighbors,
                    clearance,
                    flip_normal
                );

                if (best_projection.valid &&
                    surfaces_share_node(
                        slave_surface,
                        *surfaces[static_cast<std::size_t>(best_projection.surface_id)]
                    )) {
                    best_projection = {};
                    ++diagnostics.self_contact_rejections;
                }

                if (best_projection.valid) {
                    ++diagnostics.valid_projections;
                } else if (valid_surface_id(previous_partner->second, surfaces)) {
                    const auto& previous_surface =
                        surfaces[static_cast<std::size_t>(previous_partner->second)];

                    if (!surfaces_share_node(slave_surface, *previous_surface)) {
                        best_projection = project_on_surface(
                            previous_partner->second,
                            slave_position,
                            node_coords,
                            surfaces,
                            clearance,
                            flip_normal
                        );
                    }

                    if (best_projection.valid) {
                        ++diagnostics.valid_projections;
                    } else {
                        ++diagnostics.invalid_projections;
                    }
                }
            }

            if (!freeze_surface_partners && !best_projection.valid) {
                const auto& candidate_ids =
                    bvh.query_point(slave_position, &candidates);

                const Index candidate_count =
                    static_cast<Index>(candidate_ids.size());

                diagnostics.candidate_surfaces += candidate_count;
                diagnostics.maximum_candidates =
                    std::max(diagnostics.maximum_candidates, candidate_count);

                if (candidate_ids.empty()) {
                    ++diagnostics.zero_candidate_points;
                }

                for (ID surface_id : candidate_ids) {
                    if (!valid_surface_id(surface_id, surfaces)) {
                        continue;
                    }

                    const auto& candidate_surface =
                        surfaces[static_cast<std::size_t>(surface_id)];

                    if (surfaces_share_node(slave_surface, *candidate_surface)) {
                        ++diagnostics.self_contact_rejections;
                        continue;
                    }

                    const PatchProjection candidate_projection =
                        project_on_surface(
                            surface_id,
                            slave_position,
                            node_coords,
                            surfaces,
                            clearance,
                            flip_normal
                        );

                    if (!candidate_projection.valid) {
                        ++diagnostics.invalid_projections;
                        continue;
                    }

                    ++diagnostics.valid_projections;

                    if (candidate_projection.distance > search_radius) {
                        ++diagnostics.distance_rejections;
                        continue;
                    }

                    if (projection_is_better(candidate_projection, best_projection)) {
                        best_projection = candidate_projection;
                    }
                }
            }

            if (!best_projection.valid) {
                ++diagnostics.no_partner;

                if (previously_active || normal_multiplier > Precision(0)) {
                    ++diagnostics.active_partner_losses;
                }

                continue;
            }

            const ID best_surface_id = best_projection.surface_id;
            const Vec2 best_local = best_projection.local;
            const Precision best_distance = best_projection.distance;
            const Precision gap = best_projection.gap;

            diagnostics.maximum_closest_distance =
                std::max(diagnostics.maximum_closest_distance, best_distance);

            const auto& best_surface =
                surfaces[static_cast<std::size_t>(best_surface_id)];

            current_partners[point_id] = best_surface_id;
            current_multipliers[point_id] = normal_multiplier;
            current_gaps[point_id] = gap;
            current_characteristic_lengths[point_id] =
                std::sqrt(std::max(best_surface->area(node_coords), Precision(0)));

            const ContactLaw contact_law =
                evaluate_augmented_lagrange_law(
                    gap,
                    normal_multiplier,
                    penalty
                );

            if (!contact_law.active) {
                ++diagnostics.open_closest_partner;
                continue;
            }

            ++diagnostics.active_contacts;
            current_active_slaves.insert(point_id);

            const Precision penetration = std::max(Precision(0), -gap);
            const Precision point_force = integration_weight * contact_law.pressure;

            diagnostics.point_force_squared_sum += point_force * point_force;
            diagnostics.maximum_point_force =
                std::max(diagnostics.maximum_point_force, point_force);

            if (penetration > diagnostics.maximum_penetration) {
                diagnostics.maximum_penetration = penetration;
                diagnostics.worst_point         = point_id;
                diagnostics.worst_surface       = best_surface_id;
                diagnostics.worst_local         = best_local;
                diagnostics.worst_gap           = gap;
                diagnostics.worst_distance      = best_distance;
            }

            if (!previously_active) {
                ++diagnostics.activations;
            }

            if (previous_partner != state.partners.end() &&
                previous_partner->second != best_surface_id) {
                ++diagnostics.partner_switches;
            }

            hash_value(
                diagnostics.signature,
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(point_id)
                )
            );

            hash_value(
                diagnostics.signature,
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(best_surface_id)
                )
            );

            const ProjectionInfo projection_info =
                classify_projection(best_surface->n_nodes, best_local);

            dispatch_surface_contact_pair<NS>(
                slave_surface,
                best_surface,
                slave_local,
                best_local,
                projection_info,
                penalty,
                clearance,
                normal_multiplier,
                integration_weight,
                flip_normal,
                system_nodal_dofs,
                node_coords,
                nodal_forces,
                triplets,
                diagnostics
            );
        }
    };

    for (ID slave_surface_id : slave_surface_ids) {
        if (!valid_surface_id(slave_surface_id, surfaces)) {
            continue;
        }

        const auto& slave_surface =
            surfaces[static_cast<std::size_t>(slave_surface_id)];

        switch (slave_surface->n_nodes) {
            case 3: {
                const auto* typed =
                    dynamic_cast<const model::Surface<3>*>(slave_surface.get());
                logging::error(typed != nullptr,
                    "CONTACT: failed to cast three-node slave surface");
                process_slave_surface(*typed, slave_surface_id);
                break;
            }

            case 4: {
                const auto* typed =
                    dynamic_cast<const model::Surface<4>*>(slave_surface.get());
                logging::error(typed != nullptr,
                    "CONTACT: failed to cast four-node slave surface");
                process_slave_surface(*typed, slave_surface_id);
                break;
            }

            case 6: {
                const auto* typed =
                    dynamic_cast<const model::Surface<6>*>(slave_surface.get());
                logging::error(typed != nullptr,
                    "CONTACT: failed to cast six-node slave surface");
                process_slave_surface(*typed, slave_surface_id);
                break;
            }

            case 8: {
                const auto* typed =
                    dynamic_cast<const model::Surface<8>*>(slave_surface.get());
                logging::error(typed != nullptr,
                    "CONTACT: failed to cast eight-node slave surface");
                process_slave_surface(*typed, slave_surface_id);
                break;
            }

            default:
                logging::error(
                    false,
                    "CONTACT: unsupported slave surface with ",
                    slave_surface->n_nodes,
                    " nodes"
                );
        }
    }

    for (ID point_id : state.active_slaves) {
        if (current_active_slaves.find(point_id) == current_active_slaves.end()) {
            ++diagnostics.deactivations;
        }
    }

    const bool signature_changed =
        state.previous_signature != 0 &&
        state.previous_signature != diagnostics.signature;

    state.partners                = std::move(current_partners);
    state.active_slaves           = std::move(current_active_slaves);
    state.normal_multipliers      = std::move(current_multipliers);
    state.gaps                    = std::move(current_gaps);
    state.characteristic_lengths = std::move(current_characteristic_lengths);
    state.previous_signature      = diagnostics.signature;
    state.previous_active         = diagnostics.active_contacts;
    state.last_signature_changed  = signature_changed;

    if (!freeze_surface_partners &&
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_after_update) {
        runtime_state.trials.back().freeze_surface_partners = true;
        runtime_state.trials.back().freeze_after_update     = false;
    }

    const Precision contact_force_norm =
        std::sqrt(diagnostics.point_force_squared_sum);

    const auto assembly_end = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(assembly_end - assembly_start).count();

    const std::size_t added_triplets = triplets.size() - initial_triplet_count;

    const Precision average_candidates =
        diagnostics.slave_points > 0
            ? static_cast<Precision>(diagnostics.candidate_surfaces) /
              static_cast<Precision>(diagnostics.slave_points)
            : Precision(0);

    if (print_surface_contact_summary) {
        const auto previous_flags = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT]"
            << " call="        << runtime_state.call
            << " mode=SURFACE_IP"
            << " depth="       << runtime_state.trials.size()
            << " frozen="      << (freeze_surface_partners ? 1 : 0)
            << " points="      << diagnostics.slave_points
            << " active="      << diagnostics.active_contacts
            << " no_partner="  << diagnostics.no_partner
            << " cand_avg="    << average_candidates
            << " cand_max="    << diagnostics.maximum_candidates
            << " switches="    << diagnostics.partner_switches
            << " max_pen="     << diagnostics.maximum_penetration
            << " force_norm="  << contact_force_norm
            << " max_ip_force=" << diagnostics.maximum_point_force
            << " consistent="  << diagnostics.consistent_tangents
            << " stabilized="  << diagnostics.stabilized_tangents
            << " asym="        << diagnostics.maximum_raw_asymmetry
            << " triplets="    << added_triplets
            << " ms="          << elapsed_ms
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }
}

} // namespace fem::constraint
