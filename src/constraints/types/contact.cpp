/**
 * @file contact.cpp
 * @brief Implements frictionless node-to-surface penalty contact.
 *
 * The master surface is converted into the same topological search
 * triangulation used by CalculiX: 3/4/6/8-node faces become 1/2/4/6 search
 * triangles, neighboring triangles are connected through common edges and each
 * edge receives a common bordering plane constructed from averaged nodal master
 * normals. A slave node starts at the nearest search-triangle centroid and walks
 * across violated bordering planes until its owning search triangle is found.
 * The search triangulation is used only for robust master-face ownership; force
 * and tangent evaluation remain on the complete finite-element master face.
 *
 * Discrete contact elements are updated by Newton tangent evaluations and frozen
 * for nested line-search residual evaluations. Every retained element stores its
 * master face and representative slave area. The positive-gap cutoff is applied
 * only while regenerating the contact elements; a frozen line search therefore
 * evaluates the same discrete contact system that produced the Newton tangent.
 *
 * Following CalculiX N2F contact, contact-element regeneration is allowed through
 * the first eight Newton tangent evaluations of one increment attempt and is
 * frozen afterwards. Continuous closest-point geometry on every retained master
 * face remains current even while the discrete contact elements are frozen.
 *
 * The normal contact law follows the regularized linear node-to-face law used by
 * CalculiX. Around zero gap the linear penalty is multiplied by an arctangent
 * transition, giving a smooth force response and a very small tensile branch on
 * the open side.
 *
 * The tangent follows the analytic node-to-face linearization used by CalculiX.
 * The selected master face and representative slave area are fixed during one
 * tangent evaluation, while the closest-point coordinates, shape functions,
 * surface normal, gap and force law are differentiated consistently from the
 * closest-point orthogonality equations. The frictionless tangent is symmetrized
 * after local assembly.
 *
 * @see Contact
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 12.08.2026
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/geometry/surface/surface_interface.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fem::constraint {
namespace {

constexpr Precision geometry_tolerance     = Precision(1e-12);
constexpr Precision tangent_tolerance      = Precision(1e-14);
constexpr Precision smoothing_length_ratio = Precision(1e-6);
constexpr Precision contact_cutoff_ratio   = Precision(1e-3);
constexpr Precision search_border_ratio    = Precision(1e-3);
constexpr Precision pi                     = Precision(3.141592653589793238462643383279502884L);
constexpr Index     contact_update_iterations = 8;

constexpr bool contact_diagnostics      = true;
constexpr bool contact_pair_diagnostics = true;

/**
 * @brief Master-face ownership obtained from the search triangulation.
 */
struct MasterProjection {
    ID                           surface_id      = ID(-1);
    ID                           search_triangle = ID(-1);
    Index                        walk_steps      = 0;
    bool                         between         = false;
    model::SurfaceInterface::Ptr surface;
    Vec2                         local           = Vec2::Zero();
    Vec3                         point           = Vec3::Zero();
    Precision                    distance_sq     = std::numeric_limits<Precision>::max();
};

/**
 * @brief Current force state of one slave-node/master-face contact pair.
 */
struct PairEvaluation {
    bool            valid        = false;
    bool            active       = false;
    Precision       gap          = Precision(0);
    Precision       force_scalar = Precision(0);
    Vec2            local        = Vec2::Zero();
    Vec3            relative     = Vec3::Zero();
    Vec3            normal       = Vec3::Zero();
    Vec3            slave_force  = Vec3::Zero();
    DynamicVector   shape;
    std::vector<ID> node_ids;
};

bool valid_surface_id(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           static_cast<bool>(surfaces[static_cast<std::size_t>(surface_id)]);
}

bool surface_contains_node(const model::SurfaceInterface& surface, ID node_id) {
    for (Index local_node = 0; local_node < surface.n_nodes; ++local_node) {
        if (surface.nodes()[local_node] == node_id) {
            return true;
        }
    }
    return false;
}

/**
 * @brief CalculiX-style topological master search graph.
 *
 * The topology reproduces `triangucont` and `trianeighbor`. Geometry is rebuilt
 * from the current nodal coordinates like `updatecontpen`: master normals are
 * evaluated at every face node, accumulated over incident faces and normalized.
 * Each triangle edge then receives the common bordering plane used for the
 * triangle walk in `gencontelem_n2f`.
 */
class MasterSearchGraph {
    struct BorderPlane {
        Vec3      normal = Vec3::Zero();
        Precision offset = Precision(0);
        bool      valid  = false;

        Precision distance(const Vec3& point) const {
            return normal.dot(point) + offset;
        }
    };

    struct Triangle {
        ID                           surface_id = ID(-1);
        model::SurfaceInterface::Ptr surface;
        std::array<ID, 3>            nodes     {ID(-1), ID(-1), ID(-1)};
        std::array<ID, 3>            neighbors {ID(-1), ID(-1), ID(-1)};
        std::array<BorderPlane, 3>   borders;
        Vec3                         center = Vec3::Zero();
        bool                         valid  = false;
    };

    struct EdgeKey {
        ID first  = ID(-1);
        ID second = ID(-1);

        bool operator==(const EdgeKey& other) const {
            return first == other.first && second == other.second;
        }
    };

    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey& edge) const {
            std::size_t seed = std::hash<ID>{}(edge.first);
            seed ^= std::hash<ID>{}(edge.second) + std::size_t(0x9e3779b9) + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    struct EdgeOwner {
        ID    triangle = ID(-1);
        Index opposite = 0;
    };

public:
    struct Location {
        ID    triangle_id = ID(-1);
        ID    surface_id  = ID(-1);
        Index steps       = 0;
        bool  between     = false;
    };

private:
    std::vector<Triangle> triangles;

    static std::vector<std::array<Index, 3>> triangle_pattern(Index n_nodes) {
        switch (n_nodes) {
            case 3:
                return {{0, 1, 2}};
            case 4:
                return {{0, 1, 3}, {1, 2, 3}};
            case 6:
                return {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};
            case 8:
                return {{0, 4, 7}, {4, 1, 5}, {6, 5, 2}, {7, 6, 3}, {7, 4, 6}, {4, 5, 6}};
            default:
                return {};
        }
    }

    static EdgeKey edge_key(ID first, ID second) {
        if (second < first) {
            std::swap(first, second);
        }
        return {first, second};
    }

    static std::pair<Index, Index> edge_nodes(Index opposite) {
        switch (opposite) {
            case 0: return {1, 2};
            case 1: return {2, 0};
            default: return {0, 1};
        }
    }

    static BorderPlane build_border_plane(
        const Vec3& opposite_point,
        const Vec3& first_point,
        const Vec3& second_point,
        const Vec3& first_normal,
        const Vec3& second_normal
    ) {
        BorderPlane result;

        const Vec3 edge = second_point - first_point;
        const Precision edge_length_sq = edge.squaredNorm();
        if (!(edge_length_sq > geometry_tolerance * geometry_tolerance)) {
            return result;
        }

        Vec3 projected_first  = first_normal  - edge * (first_normal.dot(edge)  / edge_length_sq);
        Vec3 projected_second = second_normal - edge * (second_normal.dot(edge) / edge_length_sq);

        const Precision first_length  = projected_first.norm();
        const Precision second_length = projected_second.norm();
        if (!(first_length > geometry_tolerance) || !(second_length > geometry_tolerance)) {
            return result;
        }

        projected_first  /= first_length;
        projected_second /= second_length;

        Vec3 representative_normal = projected_first + projected_second;
        const Precision representative_length = representative_normal.norm();
        if (!(representative_length > geometry_tolerance)) {
            return result;
        }
        representative_normal /= representative_length;

        Vec3 plane_normal = edge.cross(representative_normal);
        const Precision plane_length = plane_normal.norm();
        if (!(plane_length > geometry_tolerance)) {
            return result;
        }
        plane_normal /= plane_length;

        Precision offset = -plane_normal.dot(first_point);

        // `straighteq3dpen` orients every bordering-plane normal outwards.
        if (plane_normal.dot(opposite_point) + offset > Precision(0)) {
            plane_normal = -plane_normal;
            offset       = -offset;
        }

        result.normal = plane_normal;
        result.offset = offset;
        result.valid  = true;
        return result;
    }

    static std::unordered_map<ID, Vec3> build_nodal_normals(
        const model::SurfaceRegion&                       master_region,
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                               node_coords
    ) {
        std::unordered_map<ID, Vec3> normals;

        for (ID surface_id : master_region) {
            if (!valid_surface_id(surface_id, surfaces)) {
                continue;
            }

            const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
            const DynamicMatrix natural = surface->node_coords_natural();
            if (natural.rows() != surface->n_nodes || natural.cols() != 2 || !natural.allFinite()) {
                continue;
            }

            for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
                const Vec2 local = natural.row(local_node).transpose();
                Vec3 normal = surface->normal(node_coords, local);
                const Precision length = normal.norm();
                if (!normal.allFinite() || !(length > geometry_tolerance)) {
                    continue;
                }

                const ID node = surface->nodes()[local_node];
                auto [entry, inserted] = normals.try_emplace(node, Vec3::Zero());
                (void)inserted;
                entry->second += normal / length;
            }
        }

        for (auto& [node, normal] : normals) {
            (void)node;
            const Precision length = normal.norm();
            if (length > geometry_tolerance) {
                normal /= length;
            }
        }

        return normals;
    }

    void build_triangles(
        const model::SurfaceRegion&                       master_region,
        const std::vector<model::SurfaceInterface::Ptr>& surfaces
    ) {
        for (ID surface_id : master_region) {
            if (!valid_surface_id(surface_id, surfaces)) {
                continue;
            }

            const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
            const auto pattern = triangle_pattern(surface->n_nodes);
            logging::error(!pattern.empty(),
                "CONTACT: master search supports only 3/4/6/8-node surfaces");

            for (const auto& local_nodes : pattern) {
                Triangle triangle;
                triangle.surface_id = surface_id;
                triangle.surface    = surface;
                for (Index local = 0; local < 3; ++local) {
                    triangle.nodes[local] = surface->nodes()[local_nodes[local]];
                }
                triangles.push_back(std::move(triangle));
            }
        }
    }

    void build_neighbors() {
        std::unordered_map<EdgeKey, EdgeOwner, EdgeKeyHash> edges;
        edges.reserve(triangles.size() * 3);

        for (ID triangle_id = 0; triangle_id < static_cast<ID>(triangles.size()); ++triangle_id) {
            auto& triangle = triangles[static_cast<std::size_t>(triangle_id)];

            for (Index opposite = 0; opposite < 3; ++opposite) {
                const auto [first_local, second_local] = edge_nodes(opposite);
                const EdgeKey key = edge_key(triangle.nodes[first_local], triangle.nodes[second_local]);
                const auto owner = edges.find(key);

                if (owner == edges.end()) {
                    edges.emplace(key, EdgeOwner{triangle_id, opposite});
                    continue;
                }

                auto& other = triangles[static_cast<std::size_t>(owner->second.triangle)];
                logging::error(other.neighbors[owner->second.opposite] < 0,
                    "CONTACT: non-manifold master search edge shared by more than two triangles");

                triangle.neighbors[opposite] = owner->second.triangle;
                other.neighbors[owner->second.opposite] = triangle_id;
            }
        }
    }

    void update_geometry(
        const model::Field&                    node_coords,
        const std::unordered_map<ID, Vec3>& nodal_normals
    ) {
        for (auto& triangle : triangles) {
            std::array<Vec3, 3> points;
            std::array<Vec3, 3> normals;

            bool valid = true;
            for (Index local = 0; local < 3; ++local) {
                points[local] = node_coords.row_vec3(triangle.nodes[local]);
                const auto normal = nodal_normals.find(triangle.nodes[local]);
                if (!points[local].allFinite() || normal == nodal_normals.end() ||
                    !normal->second.allFinite() || normal->second.norm() <= geometry_tolerance) {
                    valid = false;
                    break;
                }
                normals[local] = normal->second;
            }

            if (!valid) {
                continue;
            }

            triangle.center = (points[0] + points[1] + points[2]) / Precision(3);
            triangle.valid = true;

            for (Index opposite = 0; opposite < 3; ++opposite) {
                const auto [first_local, second_local] = edge_nodes(opposite);
                triangle.borders[opposite] = build_border_plane(
                    points[opposite],
                    points[first_local],
                    points[second_local],
                    normals[first_local],
                    normals[second_local]
                );
                triangle.valid = triangle.valid && triangle.borders[opposite].valid;
            }
        }
    }

public:
    MasterSearchGraph(
        const model::SurfaceRegion&                       master_region,
        const std::vector<model::SurfaceInterface::Ptr>& surfaces,
        const model::Field&                               node_coords
    ) {
        build_triangles(master_region, surfaces);
        build_neighbors();
        update_geometry(node_coords, build_nodal_normals(master_region, surfaces, node_coords));
    }

    std::size_t size() const {
        return triangles.size();
    }

    std::optional<Location> locate(
        ID          slave_node,
        const Vec3& point,
        Precision   border_tolerance
    ) const {
        ID seed = ID(-1);
        Precision best_distance_sq = std::numeric_limits<Precision>::max();

        for (ID triangle_id = 0; triangle_id < static_cast<ID>(triangles.size()); ++triangle_id) {
            const auto& triangle = triangles[static_cast<std::size_t>(triangle_id)];
            if (!triangle.valid || !triangle.center.allFinite() || !triangle.surface ||
                surface_contains_node(*triangle.surface, slave_node)) {
                continue;
            }

            const Precision distance_sq = (point - triangle.center).squaredNorm();
            if (!std::isfinite(distance_sq)) {
                continue;
            }

            if (seed < 0 || distance_sq < best_distance_sq - geometry_tolerance ||
                (std::abs(distance_sq - best_distance_sq) <= geometry_tolerance && triangle_id < seed)) {
                seed = triangle_id;
                best_distance_sq = distance_sq;
            }
        }

        if (seed < 0) {
            return std::nullopt;
        }

        ID current  = seed;
        ID previous = ID(-1);
        std::unordered_set<ID> visited;
        visited.reserve(std::min<std::size_t>(triangles.size(), std::size_t(100)));
        visited.insert(current);

        constexpr Index maximum_walk_steps = 100;

        for (Index step = 0; step < maximum_walk_steps; ++step) {
            const auto& triangle = triangles[static_cast<std::size_t>(current)];
            bool crossed = false;

            for (Index opposite = 0; opposite < 3; ++opposite) {
                const BorderPlane& border = triangle.borders[opposite];
                if (!border.valid) {
                    return std::nullopt;
                }

                if (border.distance(point) <= border_tolerance) {
                    continue;
                }

                const ID next = triangle.neighbors[opposite];
                if (next < 0) {
                    return std::nullopt;
                }

                const auto& next_triangle = triangles[static_cast<std::size_t>(next)];
                if (!next_triangle.surface || surface_contains_node(*next_triangle.surface, slave_node)) {
                    return std::nullopt;
                }

                // CalculiX accepts the current triangle when the walk would
                // immediately return to the previous triangle.
                if (next == previous) {
                    return Location{current, triangle.surface_id, step + 1, true};
                }

                // Longer loops are treated as a circular search path.
                if (visited.find(next) != visited.end()) {
                    return std::nullopt;
                }

                visited.insert(next);
                previous = current;
                current  = next;
                crossed  = true;
                break;
            }

            if (!crossed) {
                const auto& owner = triangles[static_cast<std::size_t>(current)];
                return Location{current, owner.surface_id, step, false};
            }
        }

        return std::nullopt;
    }
};

/**
 * @brief CalculiX-like positive tributary-area weight of one slave facet node.
 */
Precision slave_area_weight(Index n_nodes, Index local_node) {
    switch (n_nodes) {
        case 3:
            return Precision(1) / Precision(3);
        case 4:
            return Precision(1) / Precision(4);
        case 6:
            return local_node < 3 ? Precision(1) / Precision(999) : Precision(332) / Precision(999);
        case 8:
            return local_node < 4 ? Precision(0.01) : Precision(0.24);
        default:
            return Precision(0);
    }
}

Precision slave_characteristic_length(Precision slave_area) {
    return std::sqrt(std::max(slave_area, Precision(0)));
}

Precision smoothing_length(Precision slave_area) {
    return smoothing_length_ratio * slave_characteristic_length(slave_area);
}

Precision contact_cutoff(Precision slave_area) {
    return contact_cutoff_ratio * slave_characteristic_length(slave_area);
}

Precision search_border_tolerance(Precision slave_area) {
    return search_border_ratio * slave_characteristic_length(slave_area);
}

/**
 * @brief CalculiX-style arctangent transition multiplying the linear penalty.
 */
Precision smooth_contact_factor(Precision gap, Precision slave_area) {
    const Precision eps = smoothing_length(slave_area);
    if (!(eps > geometry_tolerance)) {
        return gap < Precision(0) ? Precision(1) : Precision(0);
    }

    return Precision(0.5) + std::atan(-gap / eps) / pi;
}

/**
 * @brief Derivative of the regularized scalar slave force with respect to gap.
 */
Precision smooth_force_derivative(
    Precision gap,
    Precision slave_area,
    Precision penalty
) {
    const Precision eps = smoothing_length(slave_area);
    if (!(eps > geometry_tolerance)) {
        return gap < Precision(0) ? penalty * slave_area : Precision(0);
    }

    const Precision factor = smooth_contact_factor(gap, slave_area);
    const Precision ratio  = gap / eps;
    const Precision df_dg  = -Precision(1) / (pi * eps * (Precision(1) + ratio * ratio));

    return penalty * slave_area * (factor + gap * df_dg);
}

std::unordered_map<ID, Precision> build_slave_nodal_areas(
    const model::SurfaceRegion&                       slave_region,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords
) {
    std::unordered_map<ID, Precision> areas;

    for (ID surface_id : slave_region) {
        if (!valid_surface_id(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
        const Precision area = surface->area(node_coords);

        if (!std::isfinite(area) || area <= geometry_tolerance) {
            continue;
        }

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            const Precision weight = slave_area_weight(surface->n_nodes, local_node);
            if (weight > Precision(0)) {
                areas[surface->nodes()[local_node]] += weight * area;
            }
        }
    }

    return areas;
}

std::optional<MasterProjection> project_master_surface(
    ID                                                slave_node,
    ID                                                surface_id,
    const Vec3&                                       slave_position,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords
) {
    if (!valid_surface_id(surface_id, surfaces)) {
        return std::nullopt;
    }

    const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];
    if (surface_contains_node(*surface, slave_node)) {
        return std::nullopt;
    }

    const Vec2 local = surface->global_to_local(slave_position, node_coords, true);
    if (!local.allFinite() || !surface->in_bounds(local)) {
        return std::nullopt;
    }

    const Vec3 point = surface->local_to_global(local, node_coords);
    if (!point.allFinite()) {
        return std::nullopt;
    }

    const Precision distance_sq = (slave_position - point).squaredNorm();
    if (!std::isfinite(distance_sq)) {
        return std::nullopt;
    }

    MasterProjection projection;
    projection.surface_id  = surface_id;
    projection.surface     = surface;
    projection.local       = local;
    projection.point       = point;
    projection.distance_sq = distance_sq;
    return projection;
}

std::optional<MasterProjection> find_master_projection(
    ID                                                slave_node,
    Precision                                         slave_area,
    const Vec3&                                       slave_position,
    const MasterSearchGraph&                          search_graph,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords
) {
    const auto location = search_graph.locate(slave_node, slave_position, search_border_tolerance(slave_area));
    if (!location.has_value()) {
        return std::nullopt;
    }

    auto projection = project_master_surface(
        slave_node,
        location->surface_id,
        slave_position,
        surfaces,
        node_coords
    );
    if (!projection.has_value()) {
        return std::nullopt;
    }

    projection->search_triangle = location->triangle_id;
    projection->walk_steps      = location->steps;
    projection->between         = location->between;
    return projection;
}

PairEvaluation evaluate_pair(
    ID                                  slave_node,
    Precision                           slave_area,
    const model::SurfaceInterface::Ptr& master,
    model::Field&                       node_coords,
    Precision                           penalty,
    Precision                           clearance,
    bool                                flip_normal
) {
    PairEvaluation result;

    if (!master || slave_area <= Precision(0)) {
        return result;
    }

    const Vec3 slave_position = node_coords.row_vec3(slave_node);
    const Vec2 local = master->global_to_local(slave_position, node_coords, true);
    if (!local.allFinite() || !master->in_bounds(local)) {
        return result;
    }

    const Vec3 master_position = master->local_to_global(local, node_coords);
    Vec3 normal = master->normal(node_coords, local);
    if (flip_normal) {
        normal = -normal;
    }

    const DynamicVector shape = master->shape_function(local);
    if (!slave_position.allFinite() || !master_position.allFinite() ||
        !normal.allFinite() || normal.norm() <= geometry_tolerance ||
        shape.size() != master->n_nodes || !shape.allFinite()) {
        return result;
    }
    normal.normalize();

    const Vec3 relative = slave_position - master_position;
    const Precision gap = relative.dot(normal) - clearance;
    if (!std::isfinite(gap)) {
        return result;
    }

    result.valid        = true;
    result.active       = gap <= contact_cutoff(slave_area);
    result.gap          = gap;
    result.local        = local;
    result.relative     = relative;
    result.normal       = normal;
    result.shape        = shape;
    result.force_scalar = penalty * slave_area * gap * smooth_contact_factor(gap, slave_area);
    result.slave_force  = result.force_scalar * normal;

    result.node_ids.reserve(static_cast<std::size_t>(master->n_nodes + 1));
    result.node_ids.push_back(slave_node);
    for (Index local_node = 0; local_node < master->n_nodes; ++local_node) {
        result.node_ids.push_back(master->nodes()[local_node]);
    }

    return result;
}

/**
 * @brief Returns the free natural-coordinate direction of a bounded edge projection.
 *
 * Only projections with exactly one active natural boundary constraint are
 * classified as edge projections. Corner projections deliberately fall back to
 * the existing two-parameter linearization so this change isolates edge behavior.
 */
std::optional<Vec2> bounded_edge_direction(const Vec2& local, Index master_nodes) {
    constexpr Precision boundary_tolerance = Precision(100) * geometry_tolerance;

    const Precision r = local(0);
    const Precision s = local(1);

    if (master_nodes == 3 || master_nodes == 6) {
        const bool on_r_zero = std::abs(r) <= boundary_tolerance;
        const bool on_s_zero = std::abs(s) <= boundary_tolerance;
        const bool on_sum_one = std::abs(r + s - Precision(1)) <= boundary_tolerance;
        const Index active_constraints =
            Index(on_r_zero) + Index(on_s_zero) + Index(on_sum_one);

        if (active_constraints != 1) {
            return std::nullopt;
        }
        if (on_r_zero) {
            return Vec2(Precision(0), Precision(1));
        }
        if (on_s_zero) {
            return Vec2(Precision(1), Precision(0));
        }
        return Vec2(Precision(1), Precision(-1));
    }

    if (master_nodes == 4 || master_nodes == 8) {
        const bool on_r_min = std::abs(r + Precision(1)) <= boundary_tolerance;
        const bool on_r_max = std::abs(r - Precision(1)) <= boundary_tolerance;
        const bool on_s_min = std::abs(s + Precision(1)) <= boundary_tolerance;
        const bool on_s_max = std::abs(s - Precision(1)) <= boundary_tolerance;
        const Index active_constraints =
            Index(on_r_min) + Index(on_r_max) + Index(on_s_min) + Index(on_s_max);

        if (active_constraints != 1) {
            return std::nullopt;
        }
        if (on_r_min || on_r_max) {
            return Vec2(Precision(0), Precision(1));
        }
        return Vec2(Precision(1), Precision(0));
    }

    return std::nullopt;
}

/**
 * @brief Builds the CalculiX-style analytic tangent for one fixed master face.
 */
DynamicMatrix analytic_pair_tangent(
    const PairEvaluation&               pair,
    Precision                           slave_area,
    const model::SurfaceInterface::Ptr& master,
    const model::Field&                 node_coords,
    Precision                           penalty,
    bool                                flip_normal,
    Precision&                          determinant_out,
    bool&                               valid_out
) {
    const Index master_nodes = master->n_nodes;
    const Index local_nodes  = master_nodes + 1;
    const Index local_dofs   = Index(3) * local_nodes;

    determinant_out = Precision(0);
    valid_out       = false;

    DynamicMatrix tangent = DynamicMatrix::Zero(local_dofs, local_dofs);

    const DynamicMatrix dshape  = master->shape_derivative(pair.local);
    const DynamicMatrix ddshape = master->shape_second_derivative(pair.local);
    if (dshape.rows() != master_nodes || dshape.cols() != 2 ||
        ddshape.rows() != master_nodes || ddshape.cols() != 3 ||
        !dshape.allFinite() || !ddshape.allFinite()) {
        return tangent;
    }

    DynamicMatrix coordinates(master_nodes, 3);
    for (Index node = 0; node < master_nodes; ++node) {
        coordinates.row(node) = node_coords.row_vec3(master->nodes()[node]).transpose();
    }

    const Vec3 xr  = coordinates.transpose() * dshape.col(0);
    const Vec3 xs  = coordinates.transpose() * dshape.col(1);
    const Vec3 xrr = coordinates.transpose() * ddshape.col(0);
    const Vec3 xss = coordinates.transpose() * ddshape.col(1);
    const Vec3 xrs = coordinates.transpose() * ddshape.col(2);

    const Precision orientation = flip_normal ? Precision(-1) : Precision(1);
    const Vec3 normal_vector = orientation * xr.cross(xs);
    const Precision normal_length = normal_vector.norm();
    if (!(normal_length > geometry_tolerance) || !normal_vector.allFinite()) {
        return tangent;
    }
    const Vec3 normal = normal_vector / normal_length;

    const std::optional<Vec2> edge_direction = bounded_edge_direction(pair.local, master_nodes);

    Precision c11 = Precision(0);
    Precision c12 = Precision(0);
    Precision c22 = Precision(0);
    Precision edge_inverse = Precision(0);
    Vec2 edge = Vec2::Zero();
    Vec3 edge_tangent = Vec3::Zero();

    if (edge_direction.has_value()) {
        edge = *edge_direction;
        edge_tangent = xr * edge(0) + xs * edge(1);
        const Vec3 edge_second =
            xrr * edge(0) * edge(0) +
            xss * edge(1) * edge(1) +
            Precision(2) * xrs * edge(0) * edge(1);

        const Precision edge_jacobian =
            -edge_tangent.dot(edge_tangent) + pair.relative.dot(edge_second);
        const Precision edge_scale =
            edge_tangent.squaredNorm() + std::abs(pair.relative.dot(edge_second));
        determinant_out = edge_jacobian;

        if (!std::isfinite(edge_jacobian) || !std::isfinite(edge_scale) ||
            !(edge_scale > Precision(0)) ||
            std::abs(edge_jacobian) <= geometry_tolerance * edge_scale) {
            return tangent;
        }

        edge_inverse = Precision(1) / edge_jacobian;
    } else {
        const Precision a11 = -xr.dot(xr) + pair.relative.dot(xrr);
        const Precision a12 = -xr.dot(xs) + pair.relative.dot(xrs);
        const Precision a22 = -xs.dot(xs) + pair.relative.dot(xss);
        const Precision determinant = a11 * a22 - a12 * a12;
        determinant_out = determinant;

        // The closest-point system scales with the square of the local surface
        // length. Its determinant therefore scales with length^4. Compare the
        // determinant only against a quantity with the same dimensions instead of
        // an absolute unit-scale floor, so small but well-conditioned faces are not
        // falsely classified as singular.
        const Precision determinant_scale =
            a11 * a11 + Precision(2) * a12 * a12 + a22 * a22;
        if (!std::isfinite(determinant) || !std::isfinite(determinant_scale) ||
            determinant_scale <= Precision(0) ||
            std::abs(determinant) <= geometry_tolerance * determinant_scale) {
            return tangent;
        }

        c11 =  a22 / determinant;
        c12 = -a12 / determinant;
        c22 =  a11 / determinant;
    }

    const Precision dq_dg = smooth_force_derivative(pair.gap, slave_area, penalty);

    for (Index col_node = 0; col_node < local_nodes; ++col_node) {
        const bool slave_column = col_node == 0;
        const Index master_col  = slave_column ? Index(0) : col_node - 1;

        for (Dim component = 0; component < 3; ++component) {
            Vec3 direction = Vec3::Zero();
            direction(component) = Precision(1);

            Precision dr_local = Precision(0);
            Precision ds_local = Precision(0);

            if (edge_direction.has_value()) {
                Precision edge_rhs;
                if (slave_column) {
                    edge_rhs = -edge_tangent(component);
                } else {
                    const Precision dshape_edge =
                        dshape(master_col, 0) * edge(0) + dshape(master_col, 1) * edge(1);
                    edge_rhs =
                        pair.shape(master_col) * edge_tangent(component) -
                        dshape_edge * pair.relative(component);
                }

                const Precision d_edge = edge_inverse * edge_rhs;
                dr_local = edge(0) * d_edge;
                ds_local = edge(1) * d_edge;
            } else {
                Precision b1;
                Precision b2;

                if (slave_column) {
                    b1 = -xr(component);
                    b2 = -xs(component);
                } else {
                    b1 = pair.shape(master_col) * xr(component) -
                         dshape(master_col, 0) * pair.relative(component);
                    b2 = pair.shape(master_col) * xs(component) -
                         dshape(master_col, 1) * pair.relative(component);
                }

                dr_local = c11 * b1 + c12 * b2;
                ds_local = c12 * b1 + c22 * b2;
            }

            Vec3 drelative = -xr * dr_local - xs * ds_local;
            if (slave_column) {
                drelative += direction;
            } else {
                drelative -= pair.shape(master_col) * direction;
            }

            Vec3 dxr = xrr * dr_local + xrs * ds_local;
            Vec3 dxs = xrs * dr_local + xss * ds_local;
            if (!slave_column) {
                dxr += dshape(master_col, 0) * direction;
                dxs += dshape(master_col, 1) * direction;
            }

            const Vec3 dnormal_vector = orientation * (dxr.cross(xs) + xr.cross(dxs));
            const Vec3 dnormal = (dnormal_vector - normal * normal.dot(dnormal_vector)) / normal_length;

            const Precision dgap = drelative.dot(normal) + pair.relative.dot(dnormal);
            const Vec3 dslave_force = dq_dg * dgap * normal + pair.force_scalar * dnormal;
            const Index col = Index(3) * col_node + component;

            for (Dim row_component = 0; row_component < 3; ++row_component) {
                tangent(row_component, col) = dslave_force(row_component);
            }

            for (Index row_master = 0; row_master < master_nodes; ++row_master) {
                const Precision dN = dshape(row_master, 0) * dr_local + dshape(row_master, 1) * ds_local;
                const Vec3 dmaster_force = -dN * pair.slave_force - pair.shape(row_master) * dslave_force;
                const Index row_node = row_master + 1;
                for (Dim row_component = 0; row_component < 3; ++row_component) {
                    tangent(Index(3) * row_node + row_component, col) = dmaster_force(row_component);
                }
            }
        }
    }

    tangent = Precision(0.5) * (tangent + tangent.transpose()).eval();
    valid_out = tangent.allFinite();
    return tangent;
}

void add_force(model::NodeData& nodal_forces, ID node_id, const Vec3& force) {
    for (Dim component = 0; component < 3; ++component) {
        nodal_forces(node_id, component) += force(component);
    }
}

void add_tangent_entry(
    ID            row_node,
    Dim           row_component,
    ID            col_node,
    Dim           col_component,
    Precision     value,
    SystemDofIds& system_nodal_dofs,
    TripletList&  triplets
) {
    if (std::abs(value) <= tangent_tolerance) {
        return;
    }

    const int row = system_nodal_dofs(row_node, row_component);
    const int col = system_nodal_dofs(col_node, col_component);
    if (row >= 0 && col >= 0) {
        triplets.emplace_back(row, col, value);
    }
}

} // namespace

Contact::Contact(
    model::SurfaceRegion::Ptr master,
    model::SurfaceRegion::Ptr slave,
    Precision                 penalty_stiffness,
    Precision                 contact_clearance,
    bool                      flip_master_normal
)
    : master_surfaces(std::move(master)),
      slave_surfaces(std::move(slave)),
      penalty(penalty_stiffness),
      clearance(contact_clearance),
      flip_normal(flip_master_normal) {
    logging::error(master_surfaces != nullptr && master_surfaces->size() > 0,
        "CONTACT: master surface region is empty or undefined");
    logging::error(slave_surfaces != nullptr && slave_surfaces->size() > 0,
        "CONTACT: slave surface region is empty or undefined");
    logging::error(std::isfinite(penalty) && penalty > Precision(0),
        "CONTACT: PENALTY must be finite and positive");
    logging::error(std::isfinite(clearance),
        "CONTACT: CLEARANCE must be finite");
}

void Contact::begin_update_trial() {
    const bool outer_trial = trial_stack.empty();
    trial_stack.push_back({partners, allow_partner_updates, regeneration_frozen, topology_changed, update_iterations});

    if (outer_trial) {
        regeneration_frozen = false;
        update_iterations   = 0;
    }

    allow_partner_updates = !regeneration_frozen;
    topology_changed      = false;
}

void Contact::begin_frozen_trial() {
    trial_stack.push_back({partners, allow_partner_updates, regeneration_frozen, topology_changed, update_iterations});
    allow_partner_updates = false;
    topology_changed      = false;
}

void Contact::commit_trial() {
    logging::error(!trial_stack.empty(), "CONTACT: commit without active trial");
    if (trial_stack.empty()) {
        return;
    }

    TrialState parent = std::move(trial_stack.back());
    trial_stack.pop_back();

    allow_partner_updates = parent.allow_partner_updates;
    regeneration_frozen   = parent.regeneration_frozen;
    topology_changed      = parent.topology_changed;
    update_iterations     = parent.update_iterations;
}

void Contact::rollback_trial() {
    logging::error(!trial_stack.empty(), "CONTACT: rollback without active trial");
    if (trial_stack.empty()) {
        return;
    }

    TrialState parent = std::move(trial_stack.back());
    trial_stack.pop_back();

    partners              = std::move(parent.partners);
    allow_partner_updates = parent.allow_partner_updates;
    regeneration_frozen   = parent.regeneration_frozen;
    topology_changed      = parent.topology_changed;
    update_iterations     = parent.update_iterations;
}

void Contact::freeze_partner_updates() {
    regeneration_frozen = true;
    allow_partner_updates = false;
}

bool Contact::contact_topology_changed() const {
    return topology_changed;
}

void Contact::assemble(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList*      triplets
) const {
    logging::error(master_surfaces != nullptr && slave_surfaces != nullptr,
        "CONTACT: surface regions are not defined");
    logging::error(model_data.positions != nullptr,
        "CONTACT: POSITION field is not set");

    model::Field& node_coords = *model_data.positions;
    auto& surfaces = model_data.surfaces;

    const bool update_topology = allow_partner_updates;
    std::unordered_map<ID, Precision> slave_areas;

    if (update_topology) {
        slave_areas = build_slave_nodal_areas(*slave_surfaces, surfaces, node_coords);
    } else {
        slave_areas.reserve(partners.size());
        for (const auto& [slave_node, element] : partners) {
            slave_areas.emplace(slave_node, element.slave_area);
        }
    }

    std::optional<MasterSearchGraph> search_graph;
    if (update_topology) {
        search_graph.emplace(*master_surfaces, surfaces, node_coords);
    }

    std::unordered_map<ID, ContactElement> updated_partners;
    if (update_topology) {
        updated_partners.reserve(partners.size() + 8);
    }

    Index projected_count       = 0;
    Index valid_count           = 0;
    Index active_count          = 0;
    Index between_count         = 0;
    Index tangent_failure_count = 0;
    Index maximum_walk_steps    = 0;

    Precision min_gap         = std::numeric_limits<Precision>::max();
    Precision max_gap         = -std::numeric_limits<Precision>::max();
    Precision max_penetration = Precision(0);
    Precision max_force       = Precision(0);
    Vec3 slave_resultant      = Vec3::Zero();

    ID worst_slave    = ID(-1);
    ID worst_face     = ID(-1);
    ID worst_triangle = ID(-1);
    Vec2 worst_local  = Vec2::Zero();
    Precision worst_area  = Precision(0);
    Precision worst_gap   = Precision(0);
    Precision worst_force = Precision(0);

    for (const auto& [slave_node, slave_area] : slave_areas) {
        if (slave_area <= Precision(0)) {
            continue;
        }

        const Vec3 slave_position = node_coords.row_vec3(slave_node);
        std::optional<MasterProjection> projection;

        if (update_topology) {
            projection = find_master_projection(
                slave_node,
                slave_area,
                slave_position,
                *search_graph,
                surfaces,
                node_coords
            );
        } else {
            const auto partner = partners.find(slave_node);
            if (partner == partners.end()) {
                continue;
            }

            projection = project_master_surface(
                slave_node,
                partner->second.surface_id,
                slave_position,
                surfaces,
                node_coords
            );
        }

        if (!projection.has_value()) {
            continue;
        }
        ++projected_count;
        between_count += projection->between ? 1 : 0;
        maximum_walk_steps = std::max(maximum_walk_steps, projection->walk_steps);

        const PairEvaluation pair = evaluate_pair(
            slave_node,
            slave_area,
            projection->surface,
            node_coords,
            penalty,
            clearance,
            flip_normal
        );
        if (!pair.valid) {
            continue;
        }
        ++valid_count;

        min_gap = std::min(min_gap, pair.gap);
        max_gap = std::max(max_gap, pair.gap);

        // The cutoff decides which contact elements are generated. Once a
        // Newton tangent has fixed the discrete contact system, every stored
        // element remains present throughout its frozen line-search trials and
        // the smooth positive-gap branch carries its force continuously to zero.
        if (update_topology && !pair.active) {
            continue;
        }
        ++active_count;

        if (update_topology) {
            updated_partners[slave_node] = ContactElement{projection->surface_id, slave_area};
        }

        const Precision penetration = std::max(Precision(0), -pair.gap);
        const Precision force_abs   = std::abs(pair.force_scalar);
        max_penetration = std::max(max_penetration, penetration);
        max_force       = std::max(max_force, force_abs);
        slave_resultant += pair.slave_force;

        if (penetration >= std::max(Precision(0), -worst_gap)) {
            worst_slave    = slave_node;
            worst_face     = projection->surface_id;
            worst_triangle = projection->search_triangle;
            worst_local    = pair.local;
            worst_area     = slave_area;
            worst_gap      = pair.gap;
            worst_force    = pair.force_scalar;
        }

        add_force(nodal_forces, pair.node_ids[0], pair.slave_force);
        for (Index master_node = 0; master_node < projection->surface->n_nodes; ++master_node) {
            add_force(nodal_forces,
                pair.node_ids[static_cast<std::size_t>(master_node + 1)],
                -pair.shape(master_node) * pair.slave_force);
        }

        Precision tangent_determinant = Precision(0);
        bool tangent_valid = true;
        DynamicMatrix tangent;

        if (triplets != nullptr) {
            tangent = analytic_pair_tangent(
                pair,
                slave_area,
                projection->surface,
                node_coords,
                penalty,
                flip_normal,
                tangent_determinant,
                tangent_valid
            );

            if (!tangent_valid) {
                ++tangent_failure_count;
                std::cout << "[CONTACT] tangent invalid: slave=" << slave_node
                          << " face=" << projection->surface_id
                          << " tri=" << projection->search_triangle
                          << " local=(" << pair.local(0) << ", " << pair.local(1) << ")"
                          << " gap=" << pair.gap
                          << " area=" << slave_area
                          << " detH=" << tangent_determinant << '\n';
            }

            if (contact_pair_diagnostics) {
                std::cout << "[CONTACT] pair: slave=" << slave_node
                          << " face=" << projection->surface_id
                          << " tri=" << projection->search_triangle
                          << " walk=" << projection->walk_steps
                          << " between=" << (projection->between ? "yes" : "no")
                          << " local=(" << pair.local(0) << ", " << pair.local(1) << ")"
                          << " gap=" << pair.gap
                          << " area=" << slave_area
                          << " fn=" << pair.force_scalar
                          << " detH=" << tangent_determinant
                          << " |Kt|=" << tangent.norm() << '\n';
            }

            const Index local_nodes = static_cast<Index>(pair.node_ids.size());
            for (Index row_node = 0; row_node < local_nodes; ++row_node) {
                for (Index col_node = 0; col_node < local_nodes; ++col_node) {
                    for (Dim row_component = 0; row_component < 3; ++row_component) {
                        for (Dim col_component = 0; col_component < 3; ++col_component) {
                            add_tangent_entry(
                                pair.node_ids[static_cast<std::size_t>(row_node)],
                                row_component,
                                pair.node_ids[static_cast<std::size_t>(col_node)],
                                col_component,
                                tangent(Index(3) * row_node + row_component,
                                        Index(3) * col_node + col_component),
                                system_nodal_dofs,
                                *triplets
                            );
                        }
                    }
                }
            }
        } else if (contact_pair_diagnostics) {
            std::cout << "[CONTACT] pair: slave=" << slave_node
                      << " face=" << projection->surface_id
                      << " tri=" << projection->search_triangle
                      << " walk=" << projection->walk_steps
                      << " between=" << (projection->between ? "yes" : "no")
                      << " local=(" << pair.local(0) << ", " << pair.local(1) << ")"
                      << " gap=" << pair.gap
                      << " area=" << slave_area
                      << " fn=" << pair.force_scalar << '\n';
        }
    }

    if (update_topology) {
        topology_changed = updated_partners.size() != partners.size();
        if (!topology_changed) {
            for (const auto& [slave_node, element] : updated_partners) {
                const auto previous = partners.find(slave_node);
                if (previous == partners.end() || previous->second.surface_id != element.surface_id) {
                    topology_changed = true;
                    break;
                }
            }
        }
        partners = std::move(updated_partners);

        if (triplets != nullptr) {
            ++update_iterations;
            if (update_iterations >= contact_update_iterations && !regeneration_frozen) {
                regeneration_frozen = true;
                allow_partner_updates = false;
                std::cout << "[CONTACT] topology frozen after " << update_iterations
                          << " Newton update assemblies\n";
            }
        }
    }

    const bool have_valid_gap = valid_count > 0;
    if (contact_diagnostics) {
        std::cout << "[CONTACT] summary: slaves=" << slave_areas.size()
                  << " search_triangles=" << (search_graph.has_value() ? search_graph->size() : std::size_t(0))
                  << " projected=" << projected_count
                  << " valid=" << valid_count
                  << " active=" << active_count
                  << " stored=" << partners.size()
                  << " topology=" << (update_topology ? "update" : "frozen")
                  << " changed=" << (topology_changed ? "yes" : "no")
                  << " updates=" << update_iterations
                  << " between=" << between_count
                  << " max_walk=" << maximum_walk_steps
                  << " min_gap=" << (have_valid_gap ? min_gap : Precision(0))
                  << " max_gap=" << (have_valid_gap ? max_gap : Precision(0))
                  << " max_pen=" << max_penetration
                  << " max|fn|=" << max_force
                  << " slave_resultant=(" << slave_resultant(0) << ", " << slave_resultant(1) << ", " << slave_resultant(2) << ")"
                  << " tangent_failures=" << tangent_failure_count
                  << " tangent=" << (triplets != nullptr ? "yes" : "no") << '\n';
    }

    if (contact_diagnostics && active_count > 0) {
        std::cout << "[CONTACT] worst: slave=" << worst_slave
                  << " face=" << worst_face
                  << " tri=" << worst_triangle
                  << " local=(" << worst_local(0) << ", " << worst_local(1) << ")"
                  << " gap=" << worst_gap
                  << " area=" << worst_area
                  << " fn=" << worst_force << '\n';
    }
}

} // namespace fem::constraint
