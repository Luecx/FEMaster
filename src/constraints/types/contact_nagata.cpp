/**
 * @file contact_nagata.cpp
 * @brief Implements the Nagata-smoothed master geometry used by contact.
 *
 * @author Finn Eggers
 */

#include "contact_nagata.h"

#include "../../core/logging.h"

#include <Eigen/Cholesky>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_set>

namespace fem::constraint {
namespace {

constexpr Index maximum_projection_iterations = 40;
constexpr Index maximum_line_search_iterations = 20;
constexpr Precision projection_step_tolerance = Precision(1e-12);
constexpr Precision projection_gradient_tolerance = Precision(1e-12);
constexpr Precision projection_armijo_constant = Precision(1e-4);
constexpr Precision projection_minimum_step = Precision(1e-8);
constexpr Precision projection_metric_regularization = Precision(1e-14);
constexpr Precision projection_edge_tolerance = Precision(1e-7);

class DisjointSet {
    std::vector<Index> parent_;
    std::vector<Index> rank_;

public:
    explicit DisjointSet(Index size)
        : parent_(static_cast<std::size_t>(size)),
          rank_  (static_cast<std::size_t>(size), Index(0)) {
        std::iota(parent_.begin(), parent_.end(), Index(0));
    }

    Index find(Index value) {
        Index& parent = parent_[static_cast<std::size_t>(value)];

        if (parent != value) {
            parent = find(parent);
        }

        return parent;
    }

    void unite(Index first, Index second) {
        first  = find(first);
        second = find(second);

        if (first == second) {
            return;
        }

        if (rank_[static_cast<std::size_t>(first)] <
            rank_[static_cast<std::size_t>(second)]) {
            std::swap(first, second);
        }

        parent_[static_cast<std::size_t>(second)] = first;

        if (rank_[static_cast<std::size_t>(first)] ==
            rank_[static_cast<std::size_t>(second)]) {
            ++rank_[static_cast<std::size_t>(first)];
        }
    }
};

bool valid_parent_surface(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           surfaces[static_cast<std::size_t>(surface_id)] != nullptr;
}

Precision clamp_unit(Precision value) {
    return std::max(Precision(-1), std::min(Precision(1), value));
}

} // namespace

std::size_t NagataContactSurface::EdgeKeyHash::operator()(const EdgeKey& key) const {
    std::size_t seed = std::hash<ID>{}(key.first);

    seed ^= std::hash<ID>{}(key.second) +
            std::size_t(0x9e3779b9) +
            (seed << 6) +
            (seed >> 2);

    return seed;
}

NagataContactSurface::NagataContactSurface(
    const model::SurfaceRegion&                       region,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               reference_positions,
    Precision                                         feature_angle
)
    : feature_angle_(feature_angle) {
    logging::error(
        feature_angle_ >= Precision(0) && feature_angle_ <= Precision(3.14159265358979323846),
        "CONTACT NAGATA: feature angle must lie between zero and pi"
    );

    build_patches(region, surfaces, reference_positions);
    build_neighbors();
    classify_edges();
    build_normal_groups();
    update(surfaces, reference_positions, Precision(0));
}

void NagataContactSurface::build_patches(
    const model::SurfaceRegion&                       region,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               reference_positions
) {
    Index patch_count = 0;

    for (ID surface_id : region) {
        if (!valid_parent_surface(surface_id, surfaces)) {
            continue;
        }

        switch (surfaces[static_cast<std::size_t>(surface_id)]->n_nodes) {
            case 3: patch_count += 1; break;
            case 4: patch_count += 2; break;
            case 6: patch_count += 4; break;
            case 8: patch_count += 6; break;
            default:
                logging::error(
                    false,
                    "CONTACT NAGATA: unsupported master surface with ",
                    surfaces[static_cast<std::size_t>(surface_id)]->n_nodes,
                    " nodes"
                );
        }
    }

    patches_.reserve(static_cast<std::size_t>(patch_count));

    for (ID surface_id : region) {
        if (!valid_parent_surface(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = *surfaces[static_cast<std::size_t>(surface_id)];
        const DynamicMatrix parent_local = surface.node_coords_natural();

        switch (surface.n_nodes) {
            case 3:
                add_parent_patch(surface_id, surface, parent_local, {0, 1, 2});
                break;

            case 4:
                add_parent_patch(surface_id, surface, parent_local, {0, 1, 2});
                add_parent_patch(surface_id, surface, parent_local, {0, 2, 3});
                break;

            case 6:
                add_parent_patch(surface_id, surface, parent_local, {0, 3, 5});
                add_parent_patch(surface_id, surface, parent_local, {3, 1, 4});
                add_parent_patch(surface_id, surface, parent_local, {5, 4, 2});
                add_parent_patch(surface_id, surface, parent_local, {3, 4, 5});
                break;

            case 8:
                // Four corner triangles and two triangles in the central diamond.
                add_parent_patch(surface_id, surface, parent_local, {0, 4, 7});
                add_parent_patch(surface_id, surface, parent_local, {4, 1, 5});
                add_parent_patch(surface_id, surface, parent_local, {5, 2, 6});
                add_parent_patch(surface_id, surface, parent_local, {6, 3, 7});
                add_parent_patch(surface_id, surface, parent_local, {4, 5, 6});
                add_parent_patch(surface_id, surface, parent_local, {4, 6, 7});
                break;

            default:
                break;
        }
    }

    logging::error(
        static_cast<Index>(patches_.size()) == patch_count,
        "CONTACT NAGATA: inconsistent patch count after tessellation"
    );

    update_positions(reference_positions);
}

void NagataContactSurface::add_parent_patch(
    ID                              parent_surface,
    const model::SurfaceInterface& surface,
    const DynamicMatrix&            parent_local,
    const std::array<Index, 3>&     local_nodes
) {
    NagataPatch patch;

    patch.id             = static_cast<ID>(patches_.size());
    patch.parent_surface = parent_surface;

    for (Index corner = 0; corner < 3; ++corner) {
        const Index parent_node = local_nodes[static_cast<std::size_t>(corner)];

        logging::error(
            parent_node >= 0 && parent_node < surface.n_nodes,
            "CONTACT NAGATA: invalid local node in parent-surface tessellation"
        );

        patch.nodes[static_cast<std::size_t>(corner)] = surface.nodes()[parent_node];
        patch.parent_local[static_cast<std::size_t>(corner)] =
            parent_local.row(parent_node).transpose();
    }

    patches_.push_back(std::move(patch));
}

void NagataContactSurface::build_neighbors() {
    std::unordered_map<EdgeKey, EdgeOwner, EdgeKeyHash> owners;
    owners.reserve(3 * patches_.size());

    for (NagataPatch& patch : patches_) {
        for (Index edge = 0; edge < 3; ++edge) {
            const auto corners = edge_corners(edge);

            ID first  = patch.nodes[static_cast<std::size_t>(corners[0])];
            ID second = patch.nodes[static_cast<std::size_t>(corners[1])];

            if (second < first) {
                std::swap(first, second);
            }

            const EdgeKey key{first, second};
            const auto owner = owners.find(key);

            if (owner == owners.end()) {
                owners.emplace(key, EdgeOwner{patch.id, edge});
                continue;
            }

            NagataPatch& neighbor = patches_[static_cast<std::size_t>(owner->second.patch)];

            logging::error(
                neighbor.neighbors[static_cast<std::size_t>(owner->second.edge)] == nullptr,
                "CONTACT NAGATA: non-manifold edge shared by more than two patches"
            );

            patch.neighbors[static_cast<std::size_t>(edge)] = &neighbor;
            neighbor.neighbors[static_cast<std::size_t>(owner->second.edge)] = &patch;
        }
    }
}

void NagataContactSurface::classify_edges() {
    for (NagataPatch& patch : patches_) {
        const Vec3 patch_normal = linear_normal(patch);

        for (Index edge = 0; edge < 3; ++edge) {
            NagataPatch* neighbor = patch.neighbors[static_cast<std::size_t>(edge)];

            if (neighbor == nullptr) {
                patch.edge_types[static_cast<std::size_t>(edge)] = NagataEdgeType::Boundary;
                continue;
            }

            if (patch.id > neighbor->id) {
                continue;
            }

            const auto corners = edge_corners(edge);
            const ID node_0 = patch.nodes[static_cast<std::size_t>(corners[0])];
            const ID node_1 = patch.nodes[static_cast<std::size_t>(corners[1])];
            const Index neighbor_edge = find_edge(*neighbor, node_0, node_1);

            logging::error(
                neighbor_edge >= 0,
                "CONTACT NAGATA: failed to locate reciprocal neighbor edge"
            );

            NagataEdgeType type = NagataEdgeType::Internal;

            if (patch.parent_surface != neighbor->parent_surface) {
                const Vec3 neighbor_normal = linear_normal(*neighbor);
                const Precision angle = std::acos(clamp_unit(patch_normal.dot(neighbor_normal)));

                type = angle > feature_angle_
                    ? NagataEdgeType::Sharp
                    : NagataEdgeType::Smooth;
            }

            patch.edge_types[static_cast<std::size_t>(edge)] = type;
            neighbor->edge_types[static_cast<std::size_t>(neighbor_edge)] = type;
        }
    }
}

void NagataContactSurface::build_normal_groups() {
    const Index corner_count = 3 * static_cast<Index>(patches_.size());
    DisjointSet groups(corner_count);

    for (NagataPatch& patch : patches_) {
        for (Index edge = 0; edge < 3; ++edge) {
            NagataPatch* neighbor = patch.neighbors[static_cast<std::size_t>(edge)];
            const NagataEdgeType type = patch.edge_types[static_cast<std::size_t>(edge)];

            if (neighbor == nullptr || patch.id > neighbor->id ||
                (type != NagataEdgeType::Internal && type != NagataEdgeType::Smooth)) {
                continue;
            }

            const auto corners = edge_corners(edge);

            for (Index endpoint = 0; endpoint < 2; ++endpoint) {
                const Index patch_corner = corners[static_cast<std::size_t>(endpoint)];
                const ID node = patch.nodes[static_cast<std::size_t>(patch_corner)];

                Index neighbor_corner = Index(-1);

                for (Index candidate = 0; candidate < 3; ++candidate) {
                    if (neighbor->nodes[static_cast<std::size_t>(candidate)] == node) {
                        neighbor_corner = candidate;
                        break;
                    }
                }

                logging::error(
                    neighbor_corner >= 0,
                    "CONTACT NAGATA: neighboring patches do not share the expected edge node"
                );

                groups.unite(
                    3 * static_cast<Index>(patch.id) + patch_corner,
                    3 * static_cast<Index>(neighbor->id) + neighbor_corner
                );
            }
        }
    }

    std::unordered_map<Index, Index> root_to_group;
    root_to_group.reserve(static_cast<std::size_t>(corner_count));

    for (NagataPatch& patch : patches_) {
        for (Index corner = 0; corner < 3; ++corner) {
            const Index corner_id = 3 * static_cast<Index>(patch.id) + corner;
            const Index root = groups.find(corner_id);

            const auto inserted = root_to_group.emplace(
                root,
                static_cast<Index>(normal_groups_.size())
            );

            if (inserted.second) {
                normal_groups_.emplace_back();
            }

            const Index group = inserted.first->second;
            patch.normal_groups[static_cast<std::size_t>(corner)] = group;
            normal_groups_[static_cast<std::size_t>(group)].push_back({patch.id, corner});
        }
    }
}

void NagataContactSurface::update(
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                              positions,
    Precision                                        search_radius
) {
    update_positions(positions);
    update_normals(surfaces, positions);
    update_patch_geometry();
    update_bvh(search_radius);
}

void NagataContactSurface::update_positions(const model::Field& positions) {
    for (NagataPatch& patch : patches_) {
        for (Index corner = 0; corner < 3; ++corner) {
            const ID node = patch.nodes[static_cast<std::size_t>(corner)];
            patch.positions[static_cast<std::size_t>(corner)] =
                positions.row_vec3(static_cast<Index>(node));
        }
    }
}

void NagataContactSurface::update_normals(
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                              positions
) {
    std::vector<Vec3> group_normals(normal_groups_.size(), Vec3::Zero());

    for (Index group = 0; group < static_cast<Index>(normal_groups_.size()); ++group) {
        Vec3 normal_sum = Vec3::Zero();
        std::unordered_set<ID> visited_surfaces;

        for (const CornerRef& ref : normal_groups_[static_cast<std::size_t>(group)]) {
            const NagataPatch& patch = patches_[static_cast<std::size_t>(ref.patch)];

            if (!visited_surfaces.insert(patch.parent_surface).second ||
                !valid_parent_surface(patch.parent_surface, surfaces)) {
                continue;
            }

            const auto& surface = surfaces[static_cast<std::size_t>(patch.parent_surface)];
            const Vec3 normal = surface->normal(
                positions,
                patch.parent_local[static_cast<std::size_t>(ref.corner)]
            );

            if (normal.allFinite() && normal.squaredNorm() > Precision(1e-20)) {
                normal_sum += normal.normalized();
            }
        }

        if (!normal_sum.allFinite() || normal_sum.squaredNorm() <= Precision(1e-20)) {
            const CornerRef fallback = normal_groups_[static_cast<std::size_t>(group)].front();
            normal_sum = linear_normal(patches_[static_cast<std::size_t>(fallback.patch)]);
        }

        group_normals[static_cast<std::size_t>(group)] = normal_sum.normalized();
    }

    for (NagataPatch& patch : patches_) {
        for (Index corner = 0; corner < 3; ++corner) {
            const Index group = patch.normal_groups[static_cast<std::size_t>(corner)];
            patch.normals[static_cast<std::size_t>(corner)] =
                group_normals[static_cast<std::size_t>(group)];
        }
    }
}

void NagataContactSurface::update_patch_geometry() {
    static const std::array<Vec2, 3> area_points{
        Vec2(Precision(1.0 / 6.0), Precision(1.0 / 6.0)),
        Vec2(Precision(2.0 / 3.0), Precision(1.0 / 6.0)),
        Vec2(Precision(1.0 / 6.0), Precision(2.0 / 3.0))
    };

    for (NagataPatch& patch : patches_) {
        for (Index edge = 0; edge < 3; ++edge) {
            const auto corners = edge_corners(edge);
            const Index corner_0 = corners[0];
            const Index corner_1 = corners[1];

            if (patch.edge_types[static_cast<std::size_t>(edge)] == NagataEdgeType::Sharp) {
                patch.curvature_maps[static_cast<std::size_t>(edge)].setZero();
                patch.curvatures[static_cast<std::size_t>(edge)].setZero();
                continue;
            }

            patch.curvature_maps[static_cast<std::size_t>(edge)] = curvature_map(
                patch.normals[static_cast<std::size_t>(corner_0)],
                patch.normals[static_cast<std::size_t>(corner_1)]
            );

            patch.curvatures[static_cast<std::size_t>(edge)] =
                patch.curvature_maps[static_cast<std::size_t>(edge)] *
                (
                    patch.positions[static_cast<std::size_t>(corner_1)] -
                    patch.positions[static_cast<std::size_t>(corner_0)]
                );
        }

        patch.area = Precision(0);

        for (const Vec2& local : area_points) {
            patch.area += Precision(1.0 / 6.0) * evaluate(patch, local).jacobian;
        }

        patch.bounds = BvhAabb::Aabb::invalid();

        for (const Vec3& position : patch.positions) {
            patch.bounds.expand_point(position);
        }

        for (Index edge = 0; edge < 3; ++edge) {
            const auto corners = edge_corners(edge);
            const Vec3 control = Precision(0.5) * (
                patch.positions[static_cast<std::size_t>(corners[0])] +
                patch.positions[static_cast<std::size_t>(corners[1])] -
                patch.curvatures[static_cast<std::size_t>(edge)]
            );

            patch.bounds.expand_point(control);
        }
    }
}

void NagataContactSurface::update_bvh(Precision search_radius) {
    bvh_ = BvhAabb(search_radius);

    for (const NagataPatch& patch : patches_) {
        bvh_.add_aabb(patch.id, patch.bounds);
    }

    bvh_.finalize();
}

bool NagataContactSurface::valid_patch(ID patch_id) const {
    return patch_id >= 0 && static_cast<std::size_t>(patch_id) < patches_.size();
}

NagataPatch& NagataContactSurface::patch(ID patch_id) {
    logging::error(valid_patch(patch_id), "CONTACT NAGATA: invalid patch ID");
    return patches_[static_cast<std::size_t>(patch_id)];
}

const NagataPatch& NagataContactSurface::patch(ID patch_id) const {
    logging::error(valid_patch(patch_id), "CONTACT NAGATA: invalid patch ID");
    return patches_[static_cast<std::size_t>(patch_id)];
}

const std::vector<NagataPatch>& NagataContactSurface::patches() const {
    return patches_;
}

Index NagataContactSurface::patch_count() const {
    return static_cast<Index>(patches_.size());
}

bool NagataContactSurface::valid() const {
    return bvh_.valid();
}

NagataContactSurface::Evaluation NagataContactSurface::evaluate(
    const NagataPatch& patch,
    const Vec2&        local
) const {
    const Precision r  = local(0);
    const Precision s  = local(1);
    const Precision l0 = Precision(1) - r - s;
    const Precision l1 = r;
    const Precision l2 = s;

    const Vec3& x0  = patch.positions[0];
    const Vec3& x1  = patch.positions[1];
    const Vec3& x2  = patch.positions[2];
    const Vec3& c01 = patch.curvatures[0];
    const Vec3& c12 = patch.curvatures[1];
    const Vec3& c20 = patch.curvatures[2];

    const Mat3& h01 = patch.curvature_maps[0];
    const Mat3& h12 = patch.curvature_maps[1];
    const Mat3& h20 = patch.curvature_maps[2];

    Evaluation out;

    out.position =
          l0 * x0
        + l1 * x1
        + l2 * x2
        - l0 * l1 * c01
        - l1 * l2 * c12
        - l2 * l0 * c20;

    out.first.col(0) =
          x1 - x0
        - (Precision(1) - Precision(2) * r - s) * c01
        - s * c12
        + s * c20;

    out.first.col(1) =
          x2 - x0
        + r * c01
        - r * c12
        - (Precision(1) - r - Precision(2) * s) * c20;

    out.second.col(0) = Precision(2) * c01;
    out.second.col(1) = Precision(2) * c20;
    out.second.col(2) = c01 - c12 + c20;

    const Vec3 normal_unnormalized = out.first.col(0).cross(out.first.col(1));
    out.jacobian = normal_unnormalized.norm();

    out.valid =
        out.position.allFinite() &&
        out.first.allFinite() &&
        out.second.allFinite() &&
        std::isfinite(out.jacobian) &&
        out.jacobian > Precision(1e-14);

    if (out.valid) {
        out.normal = normal_unnormalized / out.jacobian;
    }

    out.position_derivative[0] =
          l0 * Mat3::Identity()
        + l0 * l1 * h01
        - l0 * l2 * h20;

    out.position_derivative[1] =
          l1 * Mat3::Identity()
        - l0 * l1 * h01
        + l1 * l2 * h12;

    out.position_derivative[2] =
          l2 * Mat3::Identity()
        - l1 * l2 * h12
        + l0 * l2 * h20;

    out.position_derivative_r[0] =
        -Mat3::Identity() +
        (l0 - l1) * h01 +
        l2 * h20;

    out.position_derivative_r[1] =
         Mat3::Identity() -
         (l0 - l1) * h01 +
         l2 * h12;

    out.position_derivative_r[2] =
        -l2 * h12 - l2 * h20;

    out.position_derivative_s[0] =
        -Mat3::Identity() -
        l1 * h01 -
        (l0 - l2) * h20;

    out.position_derivative_s[1] =
        l1 * h01 + l1 * h12;

    out.position_derivative_s[2] =
         Mat3::Identity() -
         l1 * h12 +
         (l0 - l2) * h20;

    return out;
}

NagataContactSurface::Projection NagataContactSurface::project_on_patch(
    ID          patch_id,
    const Vec3& point,
    bool        clip
) const {
    if (!valid_patch(patch_id)) {
        return {};
    }

    const NagataPatch& target = patches_[static_cast<std::size_t>(patch_id)];
    Projection best = project_stationary(target, point, clip);

    if (!clip) {
        return best;
    }

    for (Index edge = 0; edge < 3; ++edge) {
        const Projection candidate = project_edge(target, edge, point);

        if (candidate.valid &&
            (!best.valid || candidate.distance < best.distance)) {
            best = candidate;
        }
    }

    return best;
}

NagataContactSurface::Projection NagataContactSurface::project_stationary(
    const NagataPatch& patch,
    const Vec3&        point,
    bool               require_inside
) const {
    const std::array<Vec2, 4> initial_guesses{
        Vec2(Precision(1.0 / 3.0), Precision(1.0 / 3.0)),
        Vec2(Precision(0.10), Precision(0.10)),
        Vec2(Precision(0.80), Precision(0.10)),
        Vec2(Precision(0.10), Precision(0.80))
    };

    Projection best;

    for (const Vec2& initial_guess : initial_guesses) {
        Vec2 local = initial_guess;
        Evaluation evaluation = evaluate(patch, local);

        if (!evaluation.valid) {
            continue;
        }

        Precision objective = Precision(0.5) * (evaluation.position - point).squaredNorm();

        for (Index iteration = 0; iteration < maximum_projection_iterations; ++iteration) {
            const Vec3 difference = evaluation.position - point;

            Vec2 gradient;
            gradient << difference.dot(evaluation.first.col(0)),
                        difference.dot(evaluation.first.col(1));

            if (!gradient.allFinite() || gradient.norm() < projection_gradient_tolerance) {
                break;
            }

            StaticMatrix<2, 2> hessian;
            hessian <<
                evaluation.first.col(0).dot(evaluation.first.col(0)) +
                    difference.dot(evaluation.second.col(0)),
                evaluation.first.col(0).dot(evaluation.first.col(1)) +
                    difference.dot(evaluation.second.col(2)),
                evaluation.first.col(0).dot(evaluation.first.col(1)) +
                    difference.dot(evaluation.second.col(2)),
                evaluation.first.col(1).dot(evaluation.first.col(1)) +
                    difference.dot(evaluation.second.col(1));

            Vec2 direction = Vec2::Zero();

            if (hessian.allFinite()) {
                const Eigen::LDLT<StaticMatrix<2, 2>> solver(hessian);

                if (solver.info() == Eigen::Success) {
                    direction = solver.solve(-gradient);
                }
            }

            if (!direction.allFinite() || gradient.dot(direction) >= Precision(0)) {
                StaticMatrix<2, 2> metric = evaluation.first.transpose() * evaluation.first;
                metric.diagonal().array() += projection_metric_regularization;

                const Eigen::LDLT<StaticMatrix<2, 2>> solver(metric);

                if (solver.info() == Eigen::Success) {
                    direction = solver.solve(-gradient);
                }
            }

            if (!direction.allFinite() || gradient.dot(direction) >= Precision(0)) {
                direction = -gradient;
            }

            const Precision directional_derivative = gradient.dot(direction);

            if (!std::isfinite(directional_derivative) || directional_derivative >= Precision(0)) {
                break;
            }

            Precision step = Precision(1);
            bool accepted = false;
            Vec2 accepted_local = local;
            Evaluation accepted_evaluation = evaluation;
            Precision accepted_objective = objective;

            for (Index line_search = 0;
                 line_search < maximum_line_search_iterations;
                 ++line_search) {
                const Vec2 trial_local = local + step * direction;

                if (require_inside && !in_bounds(trial_local)) {
                    step *= Precision(0.5);
                    continue;
                }

                const Evaluation trial_evaluation = evaluate(patch, trial_local);

                if (!trial_evaluation.valid) {
                    step *= Precision(0.5);
                    continue;
                }

                const Precision trial_objective =
                    Precision(0.5) * (trial_evaluation.position - point).squaredNorm();

                if (std::isfinite(trial_objective) &&
                    trial_objective <= objective +
                        projection_armijo_constant * step * directional_derivative) {
                    accepted_local      = trial_local;
                    accepted_evaluation = trial_evaluation;
                    accepted_objective  = trial_objective;
                    accepted            = true;
                    break;
                }

                step *= Precision(0.5);
            }

            if (!accepted || step < projection_minimum_step) {
                break;
            }

            const Vec2 applied_step = accepted_local - local;
            local      = accepted_local;
            evaluation = accepted_evaluation;
            objective  = accepted_objective;

            if (applied_step.norm() < projection_step_tolerance) {
                break;
            }
        }

        if (require_inside && !in_bounds(local, Precision(1e-10))) {
            continue;
        }

        const Projection candidate = make_projection(patch, local, point);

        if (candidate.valid && (!best.valid || candidate.distance < best.distance)) {
            best = candidate;
        }
    }

    return best;
}

NagataContactSurface::Projection NagataContactSurface::project_edge(
    const NagataPatch& patch,
    Index              edge,
    const Vec3&        point
) const {
    const Vec2 direction = edge_direction(edge);
    const std::array<Precision, 3> initial_guesses{
        Precision(0), Precision(0.5), Precision(1)
    };

    Projection best;

    for (Precision initial : initial_guesses) {
        Precision coordinate = initial;
        Vec2 local = edge_local(edge, coordinate);
        Evaluation evaluation = evaluate(patch, local);

        if (!evaluation.valid) {
            continue;
        }

        Precision objective = Precision(0.5) * (evaluation.position - point).squaredNorm();

        for (Index iteration = 0; iteration < maximum_projection_iterations; ++iteration) {
            const Vec3 difference = evaluation.position - point;
            const Vec3 tangent = evaluation.first * direction;

            const Precision dr = direction(0);
            const Precision ds = direction(1);
            const Vec3 curvature =
                dr * dr * evaluation.second.col(0) +
                ds * ds * evaluation.second.col(1) +
                Precision(2) * dr * ds * evaluation.second.col(2);

            const Precision gradient = difference.dot(tangent);
            const Precision hessian  = tangent.squaredNorm() + difference.dot(curvature);

            if (!std::isfinite(gradient) || std::abs(gradient) < projection_gradient_tolerance) {
                break;
            }

            Precision step_direction =
                std::isfinite(hessian) && std::abs(hessian) > Precision(1e-14)
                    ? -gradient / hessian
                    : -gradient / std::max(tangent.squaredNorm(), Precision(1e-14));

            if (!std::isfinite(step_direction) || gradient * step_direction >= Precision(0)) {
                step_direction = -gradient;
            }

            Precision step = Precision(1);
            bool accepted = false;
            Precision accepted_coordinate = coordinate;
            Evaluation accepted_evaluation = evaluation;
            Precision accepted_objective = objective;

            for (Index line_search = 0;
                 line_search < maximum_line_search_iterations;
                 ++line_search) {
                const Precision trial_coordinate = std::max(
                    Precision(0),
                    std::min(Precision(1), coordinate + step * step_direction)
                );

                const Vec2 trial_local = edge_local(edge, trial_coordinate);
                const Evaluation trial_evaluation = evaluate(patch, trial_local);

                if (!trial_evaluation.valid) {
                    step *= Precision(0.5);
                    continue;
                }

                const Precision trial_objective =
                    Precision(0.5) * (trial_evaluation.position - point).squaredNorm();

                if (std::isfinite(trial_objective) && trial_objective <= objective) {
                    accepted_coordinate = trial_coordinate;
                    accepted_evaluation = trial_evaluation;
                    accepted_objective  = trial_objective;
                    accepted            = true;
                    break;
                }

                step *= Precision(0.5);
            }

            if (!accepted) {
                break;
            }

            const Precision applied_step = accepted_coordinate - coordinate;
            coordinate = accepted_coordinate;
            evaluation = accepted_evaluation;
            objective  = accepted_objective;

            if (std::abs(applied_step) < projection_step_tolerance) {
                break;
            }
        }

        Projection candidate = make_projection(
            patch,
            edge_local(edge, coordinate),
            point
        );

        if (candidate.valid && (!best.valid || candidate.distance < best.distance)) {
            best = candidate;
        }
    }

    return best;
}

NagataContactSurface::Projection NagataContactSurface::make_projection(
    const NagataPatch& patch,
    const Vec2&        local,
    const Vec3&        point
) const {
    Projection projection;
    const Evaluation evaluation = evaluate(patch, local);

    if (!evaluation.valid || !local.allFinite()) {
        return projection;
    }

    projection.valid    = true;
    projection.patch    = &patch;
    projection.local    = local;
    projection.position = evaluation.position;
    projection.normal   = evaluation.normal;
    projection.distance = (point - evaluation.position).norm();

    classify_feature(projection);
    return projection;
}

NagataContactSurface::Projection NagataContactSurface::walk(
    ID          start_patch,
    const Vec3& point
) const {
    if (!valid_patch(start_patch)) {
        return {};
    }

    const NagataPatch* current = &patches_[static_cast<std::size_t>(start_patch)];
    std::unordered_set<ID> visited;
    visited.reserve(std::min<std::size_t>(patches_.size(), std::size_t(32)));

    const Index maximum_steps = static_cast<Index>(
        std::min<std::size_t>(patches_.size(), std::size_t(64))
    );

    for (Index step = 0; step < maximum_steps; ++step) {
        if (current == nullptr || !visited.insert(current->id).second) {
            break;
        }

        const Projection projection = project_stationary(*current, point, false);

        if (!projection.valid) {
            break;
        }

        if (in_bounds(projection.local, Precision(1e-9))) {
            return projection;
        }

        const ID edge = crossed_edge(projection.local);

        if (edge < 0) {
            break;
        }

        current = current->neighbors[static_cast<std::size_t>(edge)];
    }

    return {};
}

const std::vector<ID>& NagataContactSurface::query_point(
    const Vec3& point,
    std::vector<ID>* buffer
) const {
    return bvh_.query_point(point, buffer);
}

bool NagataContactSurface::in_bounds(const Vec2& local, Precision tolerance) {
    const Precision r = local(0);
    const Precision s = local(1);

    return r >= -tolerance &&
           s >= -tolerance &&
           r + s <= Precision(1) + tolerance;
}

ID NagataContactSurface::crossed_edge(const Vec2& local) {
    const Precision r = local(0);
    const Precision s = local(1);

    Precision largest_violation = Precision(1e-10);
    ID edge = ID(-1);

    const auto update = [&](Precision violation, ID candidate) {
        if (violation > largest_violation) {
            largest_violation = violation;
            edge = candidate;
        }
    };

    update(-s, 0);
    update(r + s - Precision(1), 1);
    update(-r, 2);

    return edge;
}

Vec2 NagataContactSurface::edge_direction(Index edge) {
    switch (edge) {
        case 0: return Vec2(Precision(1), Precision(0));
        case 1: return Vec2(Precision(-1), Precision(1));
        default: return Vec2(Precision(0), Precision(-1));
    }
}

std::array<Index, 2> NagataContactSurface::edge_corners(Index edge) {
    switch (edge) {
        case 0: return {0, 1};
        case 1: return {1, 2};
        default: return {2, 0};
    }
}

Index NagataContactSurface::find_edge(
    const NagataPatch& patch,
    ID                 node_0,
    ID                 node_1
) {
    for (Index edge = 0; edge < 3; ++edge) {
        const auto corners = edge_corners(edge);
        const ID first  = patch.nodes[static_cast<std::size_t>(corners[0])];
        const ID second = patch.nodes[static_cast<std::size_t>(corners[1])];

        if ((first == node_0 && second == node_1) ||
            (first == node_1 && second == node_0)) {
            return edge;
        }
    }

    return Index(-1);
}

Vec3 NagataContactSurface::linear_normal(const NagataPatch& patch) {
    const Vec3 normal =
        (patch.positions[1] - patch.positions[0]).cross(
            patch.positions[2] - patch.positions[0]
        );

    logging::error(
        normal.allFinite() && normal.squaredNorm() > Precision(1e-20),
        "CONTACT NAGATA: degenerate tessellated master patch"
    );

    return normal.normalized();
}

Mat3 NagataContactSurface::curvature_map(
    const Vec3& normal_0,
    const Vec3& normal_1
) {
    const Vec3 n0 = normal_0.normalized();
    const Vec3 n1 = normal_1.normalized();

    StaticMatrix<2, 3> constraint;
    constraint.row(0) = n0.transpose();
    constraint.row(1) = n1.transpose();

    StaticMatrix<2, 3> right_map;
    right_map.row(0) =  n0.transpose();
    right_map.row(1) = -n1.transpose();

    Eigen::JacobiSVD<StaticMatrix<2, 3>> solver(
        constraint,
        Eigen::ComputeFullU | Eigen::ComputeFullV
    );

    solver.setThreshold(
        Precision(64) * std::numeric_limits<Precision>::epsilon()
    );

    const Mat3 map = solver.solve(right_map);
    return map.allFinite() ? map : Mat3::Zero();
}

Vec2 NagataContactSurface::edge_local(Index edge, Precision coordinate) {
    switch (edge) {
        case 0:
            return Vec2(coordinate, Precision(0));
        case 1:
            return Vec2(Precision(1) - coordinate, coordinate);
        default:
            return Vec2(Precision(0), Precision(1) - coordinate);
    }
}

void NagataContactSurface::classify_feature(Projection& projection) {
    const Precision r = projection.local(0);
    const Precision s = projection.local(1);
    const Precision t = Precision(1) - r - s;

    const bool on_edge_0 = std::abs(s) <= projection_edge_tolerance;
    const bool on_edge_1 = std::abs(t) <= projection_edge_tolerance;
    const bool on_edge_2 = std::abs(r) <= projection_edge_tolerance;

    const Index edge_count =
        (on_edge_0 ? 1 : 0) +
        (on_edge_1 ? 1 : 0) +
        (on_edge_2 ? 1 : 0);

    if (edge_count >= 2) {
        projection.feature = NagataProjectionFeature::Corner;

        if (r <= projection_edge_tolerance && s <= projection_edge_tolerance) {
            projection.feature_index = 0;
        } else if (s <= projection_edge_tolerance && t <= projection_edge_tolerance) {
            projection.feature_index = 1;
        } else {
            projection.feature_index = 2;
        }

        return;
    }

    if (on_edge_0) {
        projection.feature       = NagataProjectionFeature::Edge;
        projection.feature_index = 0;
        return;
    }

    if (on_edge_1) {
        projection.feature       = NagataProjectionFeature::Edge;
        projection.feature_index = 1;
        return;
    }

    if (on_edge_2) {
        projection.feature       = NagataProjectionFeature::Edge;
        projection.feature_index = 2;
        return;
    }

    projection.feature       = NagataProjectionFeature::Interior;
    projection.feature_index = Index(-1);
}

} // namespace fem::constraint
