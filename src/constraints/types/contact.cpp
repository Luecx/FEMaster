/**
 * @file contact.cpp
 * @brief Implements frictionless node-to-Nagata-surface augmented-Lagrange contact.
 *
 * The master surface set is tessellated into topologically connected Nagata
 * triangles. New or released slave points acquire their master patch through a
 * bounded closest-point search. Established active or augmented points are
 * tracked only across neighboring Nagata patches, so ownership changes occur
 * when the unconstrained projection crosses an actual patch edge.
 *
 * The inner Newton solve keeps the selected master patch and normal multiplier
 * fixed while reprojecting the natural coordinates on that patch. Nodal normals
 * and Nagata edge maps are updated for every assembly but treated as lagged
 * quantities within the local contact linearization.
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

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

constexpr Precision augmentation_gap_relative_tolerance = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

constexpr bool  print_contact_summary = true;
constexpr bool  print_contact_details = false;
constexpr Index maximum_detail_prints = 16;

constexpr Precision nagata_feature_angle =
    Precision(0.78539816339744830962); // 45 degrees

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
    Index self_contact_rejections = 0;
    Index active_partner_losses   = 0;
    Index retained_penetrations   = 0;
    Index consistent_tangents     = 0;
    Index stabilized_tangents     = 0;

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

    ID             worst_slave    = ID(-1);
    ID             worst_surface  = ID(-1);
    Vec2           worst_local    = Vec2::Zero();
    ProjectionMode worst_mode     = ProjectionMode::Interior;
    Precision      worst_gap      = Precision(0);
    Precision      worst_distance = Precision(0);

    std::uint64_t signature = 1469598103934665603ull;
};

struct PatchProjection {
    bool valid = false;

    ID   patch_id = ID(-1);
    Vec2 local    = Vec2::Zero();

    Precision distance = std::numeric_limits<Precision>::max();
    Precision gap      = std::numeric_limits<Precision>::max();

    NagataProjectionFeature feature = NagataProjectionFeature::Interior;
    Index feature_index = Index(-1);
};

const char* projection_mode_name(ProjectionMode mode) {
    switch (mode) {
        case ProjectionMode::Interior: return "INTERIOR";
        case ProjectionMode::Edge:     return "EDGE";
        case ProjectionMode::Corner:   return "CORNER";
    }

    return "UNKNOWN";
}

ProjectionInfo classify_projection(const PatchProjection& projection) {
    ProjectionInfo info;

    switch (projection.feature) {
        case NagataProjectionFeature::Interior:
            info.mode = ProjectionMode::Interior;
            break;

        case NagataProjectionFeature::Edge:
            info.mode = ProjectionMode::Edge;
            info.edge_direction = NagataContactSurface::edge_direction(
                projection.feature_index
            );
            break;

        case NagataProjectionFeature::Corner:
            info.mode = ProjectionMode::Corner;
            break;
    }

    return info;
}

bool valid_parent_surface(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           surfaces[static_cast<std::size_t>(surface_id)] != nullptr;
}

bool parent_surface_contains_node(
    const NagataPatch&                                patch,
    ID                                                node_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    if (!valid_parent_surface(patch.parent_surface, surfaces)) {
        return false;
    }

    const auto& surface = surfaces[static_cast<std::size_t>(patch.parent_surface)];

    for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
        if (surface->nodes()[local_node] == node_id) {
            return true;
        }
    }

    return false;
}

PatchProjection evaluate_projection(
    const NagataContactSurface::Projection& geometry_projection,
    const Vec3&                             slave_position,
    Precision                               clearance,
    bool                                    flip_normal
) {
    PatchProjection projection;

    if (!geometry_projection.valid || geometry_projection.patch == nullptr) {
        return projection;
    }

    Vec3 normal = geometry_projection.normal;

    if (flip_normal) {
        normal = -normal;
    }

    const Vec3 difference = slave_position - geometry_projection.position;
    const Precision gap = difference.dot(normal) - clearance;

    if (!normal.allFinite() || !difference.allFinite() ||
        !std::isfinite(geometry_projection.distance) || !std::isfinite(gap)) {
        return projection;
    }

    projection.valid         = true;
    projection.patch_id      = geometry_projection.patch->id;
    projection.local         = geometry_projection.local;
    projection.distance      = geometry_projection.distance;
    projection.gap           = gap;
    projection.feature       = geometry_projection.feature;
    projection.feature_index = geometry_projection.feature_index;

    return projection;
}

PatchProjection project_on_patch(
    NagataContactSurface& geometry,
    ID                    patch_id,
    const Vec3&           slave_position,
    Precision             clearance,
    bool                  flip_normal
) {
    return evaluate_projection(
        geometry.project_on_patch(patch_id, slave_position, true),
        slave_position,
        clearance,
        flip_normal
    );
}

PatchProjection walk_master_patch(
    NagataContactSurface& geometry,
    ID                    start_patch,
    const Vec3&           slave_position,
    Precision             clearance,
    bool                  flip_normal
) {
    return evaluate_projection(
        geometry.walk(start_patch, slave_position),
        slave_position,
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
           candidate.patch_id < current.patch_id;
}

void hash_value(std::uint64_t& signature, std::uint64_t value) {
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
    auto& surfaces = model_data.surfaces;

    for (ID surface_id : *slave_surfaces) {
        if (!valid_parent_surface(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

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
        if (!valid_parent_surface(surface_id, surfaces)) {
            continue;
        }

        const auto& surface = surfaces[static_cast<std::size_t>(surface_id)];

        if (surface->n_nodes <= 0) {
            continue;
        }

        const Precision surface_area = surface->area(node_coords);

        logging::error(
            std::isfinite(surface_area) && surface_area > Precision(0),
            "CONTACT: invalid slave surface area for surface ",
            surface_id
        );

        const Precision nodal_area =
            surface_area / static_cast<Precision>(surface->n_nodes);

        for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
            areas[surface->nodes()[local_node]] += nodal_area;
        }
    }

    return areas;
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

template<Index N>
void scatter_contact_tangent(
    const std::array<ID, N>&          node_ids,
    const StaticMatrix<3 * N, 3 * N>& tangent,
    SystemDofIds&                     system_nodal_dofs,
    TripletList&                      triplets
) {
    constexpr Index local_dofs = 3 * N;
    std::array<int, local_dofs> global_dofs{};

    for (Index local_node = 0; local_node < N; ++local_node) {
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

void assemble_contact_pair(
    const NagataContactSurface& geometry,
    const NagataPatch&          patch,
    ID                          slave_node_id,
    const Vec2&                 master_local,
    const ProjectionInfo&       projection_info,
    Precision                   penalty,
    Precision                   clearance,
    Precision                   normal_multiplier,
    Precision                   slave_weight,
    bool                        flip_normal,
    SystemDofIds&               system_nodal_dofs,
    const model::Field&         node_coords,
    model::NodeData&            nodal_forces,
    TripletList&                triplets,
    ContactDiagnostics&         diagnostics
) {
    constexpr Index local_nodes = 4;
    constexpr Index local_dofs  = 3 * local_nodes;

    using LocalMatrix   = StaticMatrix<local_dofs, local_dofs>;
    using LocalRow      = StaticMatrix<1, local_dofs>;
    using LocalMap      = StaticMatrix<3, local_dofs>;
    using ProjectionMap = StaticMatrix<2, local_dofs>;

    const NagataContactSurface::Evaluation evaluation =
        geometry.evaluate(patch, master_local);

    logging::error(
        evaluation.valid,
        "CONTACT: singular Nagata master-patch Jacobian"
    );

    diagnostics.minimum_surface_jacobian = std::min(
        diagnostics.minimum_surface_jacobian,
        evaluation.jacobian
    );

    const Vec3 slave_position =
        node_coords.row_vec3(static_cast<Index>(slave_node_id));

    const Vec3 master_position = evaluation.position;
    const Vec3 difference      = slave_position - master_position;

    const Vec3 tangent_r  = evaluation.first.col(0);
    const Vec3 tangent_s  = evaluation.first.col(1);
    const Vec3 tangent_rr = evaluation.second.col(0);
    const Vec3 tangent_ss = evaluation.second.col(1);
    const Vec3 tangent_rs = evaluation.second.col(2);

    const Vec3 normal_base = evaluation.normal;
    const Vec3 normal = flip_normal ? -normal_base : normal_base;
    const Precision gap = difference.dot(normal) - clearance;

    std::array<ID, local_nodes> local_node_ids{
        slave_node_id,
        patch.nodes[0],
        patch.nodes[1],
        patch.nodes[2]
    };

    const ContactLaw contact_law = evaluate_augmented_lagrange_law(
        gap,
        normal_multiplier,
        penalty
    );

    if (!contact_law.active) {
        return;
    }

    const Precision contact_scale = slave_weight * penalty;
    const Vec3 slave_force = contact_scale * contact_law.value * normal;

    add_translational_force(nodal_forces, slave_node_id, slave_force);

    for (Index master_node = 0; master_node < 3; ++master_node) {
        add_translational_force(
            nodal_forces,
            patch.nodes[static_cast<std::size_t>(master_node)],
            -evaluation.position_derivative[static_cast<std::size_t>(master_node)].transpose() *
                slave_force
        );
    }

    LocalMap direct_difference = LocalMap::Zero();
    direct_difference.template block<3, 3>(0, 0) = Mat3::Identity();

    LocalMap direct_tangent_r = LocalMap::Zero();
    LocalMap direct_tangent_s = LocalMap::Zero();

    for (Index master_node = 0; master_node < 3; ++master_node) {
        const Index column = 3 * (master_node + 1);

        direct_difference.template block<3, 3>(0, column) =
            -evaluation.position_derivative[static_cast<std::size_t>(master_node)];

        direct_tangent_r.template block<3, 3>(0, column) =
            evaluation.position_derivative_r[static_cast<std::size_t>(master_node)];

        direct_tangent_s.template block<3, 3>(0, column) =
            evaluation.position_derivative_s[static_cast<std::size_t>(master_node)];
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

        if (projection_hessian.allFinite()) {
            diagnostics.minimum_interior_hessian_ratio = std::min(
                diagnostics.minimum_interior_hessian_ratio,
                projection_hessian.determinant() / (scale * scale)
            );
        }

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
        const Vec3 edge_tangent   = evaluation.first * edge_direction;

        const Precision dr = edge_direction(0);
        const Precision ds = edge_direction(1);

        const Vec3 edge_curvature =
            dr * dr * tangent_rr +
            ds * ds * tangent_ss +
            Precision(2) * dr * ds * tangent_rs;

        const LocalMap direct_edge_tangent =
            dr * direct_tangent_r + ds * direct_tangent_s;

        const LocalRow edge_rhs =
            edge_tangent.transpose() * direct_difference +
            difference.transpose() * direct_edge_tangent;

        const Precision edge_hessian =
            edge_tangent.squaredNorm() - difference.dot(edge_curvature);

        const Precision edge_scale =
            std::max(edge_tangent.squaredNorm(), Precision(1e-16));

        if (std::isfinite(edge_hessian)) {
            diagnostics.minimum_edge_hessian_ratio = std::min(
                diagnostics.minimum_edge_hessian_ratio,
                edge_hessian / edge_scale
            );
        }

        if (std::isfinite(edge_hessian) &&
            edge_hessian > Precision(1e-12) * edge_scale) {
            projection_derivative =
                edge_direction * (edge_rhs / edge_hessian);

            projection_linearization_valid = projection_derivative.allFinite();
        }
    }

    if (projection_linearization_valid &&
        projection_derivative.norm() > maximum_projection_derivative_norm) {
        projection_derivative.setZero();
        projection_linearization_valid = false;
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

    diagnostics.maximum_projection_derivative = std::max(
        diagnostics.maximum_projection_derivative,
        projection_derivative.norm()
    );

    LocalMatrix contact_tangent = LocalMatrix::Zero();

    const bool use_consistent_tangent =
        projection_info.mode == ProjectionMode::Interior &&
        projection_linearization_valid;

    if (use_consistent_tangent) {
        const LocalMap difference_derivative =
            direct_difference - evaluation.first * projection_derivative;

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
            normal_projector / evaluation.jacobian *
            (
                -skew(tangent_s) * tangent_r_derivative +
                 skew(tangent_r) * tangent_s_derivative
            );

        const LocalMap normal_derivative =
            flip_normal ? -normal_derivative_base : normal_derivative_base;

        const LocalRow gap_derivative =
            normal.transpose() * difference_derivative +
            difference.transpose() * normal_derivative;

        const LocalMap slave_force_derivative =
            contact_scale *
            (
                normal * (contact_law.derivative * gap_derivative) +
                contact_law.value * normal_derivative
            );

        contact_tangent.template block<3, local_dofs>(0, 0) =
            slave_force_derivative;

        for (Index master_node = 0; master_node < 3; ++master_node) {
            const Mat3& position_map =
                evaluation.position_derivative[static_cast<std::size_t>(master_node)];

            const Vec3 derivative_r_force =
                evaluation.position_derivative_r[static_cast<std::size_t>(master_node)].transpose() *
                slave_force;

            const Vec3 derivative_s_force =
                evaluation.position_derivative_s[static_cast<std::size_t>(master_node)].transpose() *
                slave_force;

            const LocalMap position_map_force_derivative =
                derivative_r_force * projection_derivative.row(0) +
                derivative_s_force * projection_derivative.row(1);

            contact_tangent.template block<3, local_dofs>(
                3 * (master_node + 1),
                0
            ) =
                -position_map.transpose() * slave_force_derivative -
                position_map_force_derivative;
        }

        ++diagnostics.consistent_tangents;
    } else {
        const LocalRow gap_derivative = normal.transpose() * direct_difference;

        contact_tangent =
            contact_scale * contact_law.derivative *
            gap_derivative.transpose() * gap_derivative;

        ++diagnostics.stabilized_tangents;
    }

    logging::error(
        contact_tangent.allFinite(),
        "CONTACT: non-finite local Nagata contact tangent"
    );

    const Precision tangent_norm = contact_tangent.norm();
    const Precision raw_asymmetry =
        (contact_tangent - contact_tangent.transpose()).norm() /
        std::max(tangent_norm, std::numeric_limits<Precision>::epsilon());

    diagnostics.maximum_raw_asymmetry = std::max(
        diagnostics.maximum_raw_asymmetry,
        raw_asymmetry
    );

    diagnostics.maximum_local_tangent_norm = std::max(
        diagnostics.maximum_local_tangent_norm,
        tangent_norm
    );

    scatter_contact_tangent<local_nodes>(
        local_node_ids,
        contact_tangent,
        system_nodal_dofs,
        triplets
    );
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

void Contact::begin_trial(bool freeze_partners, bool freeze_after_update) const {
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

    AssemblyState accepted_state = std::move(runtime_state.trials.back().state);
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

bool Contact::update_augmented_lagrange() const {
    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    if (state.last_signature_changed) {
        state.last_augmentation_changed = false;

        if (print_contact_summary) {
            const auto previous_flags = std::cout.flags();
            const auto previous_precision = std::cout.precision();

            std::cout
                << std::scientific
                << std::setprecision(3)
                << "[CONTACT_AUGMENT]"
                << " deferred=1"
                << " updated=0"
                << " max_dlambda=" << Precision(0)
                << " max_pen="     << Precision(0)
                << '\n';

            std::cout.flags(previous_flags);
            std::cout.precision(previous_precision);
        }

        return false;
    }

    Index     updated_points      = 0;
    Precision maximum_penetration = Precision(0);
    Precision maximum_change      = Precision(0);

    for (const auto& [slave_node_id, patch_id] : state.partners) {
        (void) patch_id;

        const auto gap_it = state.gaps.find(slave_node_id);

        if (gap_it == state.gaps.end() || !std::isfinite(gap_it->second)) {
            continue;
        }

        const Precision gap = gap_it->second;
        const auto length_it = state.characteristic_lengths.find(slave_node_id);

        const Precision characteristic_length =
            length_it != state.characteristic_lengths.end()
                ? std::max(length_it->second, Precision(0))
                : Precision(0);

        const Precision gap_tolerance = std::max(
            augmentation_gap_absolute_tolerance,
            augmentation_gap_relative_tolerance * characteristic_length
        );

        Precision& normal_multiplier = state.normal_multipliers[slave_node_id];
        const Precision old_multiplier = std::max(normal_multiplier, Precision(0));
        Precision new_multiplier = old_multiplier;

        if (gap < -gap_tolerance ||
            (gap > gap_tolerance && old_multiplier > Precision(0))) {
            new_multiplier = std::max(
                Precision(0),
                old_multiplier - penalty * gap
            );
        }

        const Precision multiplier_change =
            std::abs(new_multiplier - old_multiplier);

        const Precision multiplier_scale = std::max({
            old_multiplier,
            new_multiplier,
            penalty * gap_tolerance,
            Precision(1)
        });

        const Precision multiplier_tolerance = std::max(
            augmentation_multiplier_absolute_tolerance,
            augmentation_multiplier_relative_tolerance * multiplier_scale
        );

        normal_multiplier = new_multiplier;
        maximum_penetration = std::max(
            maximum_penetration,
            std::max(Precision(0), -gap)
        );
        maximum_change = std::max(maximum_change, multiplier_change);

        if (multiplier_change > multiplier_tolerance) {
            ++updated_points;
        }
    }

    state.last_augmentation_changed = updated_points > 0;

    if (print_contact_summary) {
        const auto previous_flags = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT_AUGMENT]"
            << " updated="      << updated_points
            << " max_dlambda="  << maximum_change
            << " max_pen="      << maximum_penetration
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }

    return state.last_augmentation_changed;
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

    logging::error(
        master_surfaces != nullptr,
        "CONTACT: master surface region is not set"
    );

    const model::Field& node_coords = *model_data.positions;
    auto& surfaces = model_data.surfaces;

    const auto assembly_start = std::chrono::steady_clock::now();
    const std::size_t initial_triplet_count = triplets.size();

    ContactDiagnostics diagnostics;
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

    if (!runtime_state.master_geometry) {
        const model::Field& reference_positions =
            model_data.positions_reference != nullptr
                ? *model_data.positions_reference
                : node_coords;

        runtime_state.master_geometry = std::make_unique<NagataContactSurface>(
            *master_surfaces,
            surfaces,
            reference_positions,
            nagata_feature_angle
        );
    }

    NagataContactSurface& master_geometry = *runtime_state.master_geometry;
    master_geometry.update(surfaces, node_coords, search_radius);

    if (print_contact_summary && !master_geometry.valid()) {
        std::cout
            << "[CONTACT]"
            << " call=" << runtime_state.call
            << " invalid_bvh=1"
            << '\n';
    }

    if (!runtime_state.slave_nodes_initialized) {
        runtime_state.slave_node_ids = collect_slave_nodes(
            slave_nodes,
            slave_surfaces,
            model_data
        );

        runtime_state.slave_nodes_initialized = true;
    }

    const auto& slave_node_ids = runtime_state.slave_node_ids;
    diagnostics.slave_nodes = static_cast<Index>(slave_node_ids.size());

    if (!runtime_state.slave_weights_initialized) {
        runtime_state.slave_tributary_areas = build_slave_tributary_areas(
            slave_surfaces,
            model_data,
            node_coords
        );

        runtime_state.slave_weights_initialized = true;
    }

    const auto& slave_tributary_areas = runtime_state.slave_tributary_areas;

    std::unordered_map<ID, ID> current_partners;
    std::unordered_set<ID>     current_active_slaves;

    current_partners.reserve(slave_node_ids.size());
    current_active_slaves.reserve(slave_node_ids.size());

    std::vector<ID> candidates;
    candidates.reserve(64);

    Index detail_prints = 0;

    std::unordered_map<ID, Precision> current_multipliers;
    std::unordered_map<ID, Precision> current_gaps;
    std::unordered_map<ID, Precision> current_characteristic_lengths;

    current_multipliers.reserve(slave_node_ids.size());
    current_gaps.reserve(slave_node_ids.size());
    current_characteristic_lengths.reserve(slave_node_ids.size());

    for (ID slave_node_id : slave_node_ids) {
        const Vec3 slave_position =
            node_coords.row_vec3(static_cast<Index>(slave_node_id));

        const auto previous_partner = state.partners.find(slave_node_id);
        const auto previous_multiplier = state.normal_multipliers.find(slave_node_id);

        const Precision normal_multiplier =
            previous_multiplier != state.normal_multipliers.end()
                ? std::max(previous_multiplier->second, Precision(0))
                : Precision(0);

        const bool previously_active =
            state.active_slaves.find(slave_node_id) != state.active_slaves.end();

        const bool track_established_partner =
            previous_partner != state.partners.end() &&
            (previously_active || normal_multiplier > Precision(0));

        PatchProjection best_projection;

        if (freeze_surface_partners) {
            if (previous_partner != state.partners.end() &&
                master_geometry.valid_patch(previous_partner->second)) {
                const NagataPatch& fixed_patch =
                    master_geometry.patch(previous_partner->second);

                if (parent_surface_contains_node(fixed_patch, slave_node_id, surfaces)) {
                    ++diagnostics.self_contact_rejections;
                } else {
                    best_projection = project_on_patch(
                        master_geometry,
                        fixed_patch.id,
                        slave_position,
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
                master_geometry,
                previous_partner->second,
                slave_position,
                clearance,
                flip_normal
            );

            if (best_projection.valid) {
                ++diagnostics.valid_projections;
            } else if (master_geometry.valid_patch(previous_partner->second)) {
                best_projection = project_on_patch(
                    master_geometry,
                    previous_partner->second,
                    slave_position,
                    clearance,
                    flip_normal
                );

                if (best_projection.valid) {
                    ++diagnostics.valid_projections;
                    ++diagnostics.retained_penetrations;
                } else {
                    ++diagnostics.invalid_projections;
                }
            }
        }

        if (!freeze_surface_partners && !best_projection.valid) {
            const auto& candidate_ids =
                master_geometry.query_point(slave_position, &candidates);

            const Index candidate_count = static_cast<Index>(candidate_ids.size());
            diagnostics.candidate_surfaces += candidate_count;
            diagnostics.maximum_candidates = std::max(
                diagnostics.maximum_candidates,
                candidate_count
            );

            if (candidate_ids.empty()) {
                ++diagnostics.zero_candidate_slaves;
            }

            for (ID patch_id : candidate_ids) {
                if (!master_geometry.valid_patch(patch_id)) {
                    continue;
                }

                const NagataPatch& candidate_patch = master_geometry.patch(patch_id);

                if (parent_surface_contains_node(candidate_patch, slave_node_id, surfaces)) {
                    ++diagnostics.self_contact_rejections;
                    continue;
                }

                const PatchProjection candidate_projection = project_on_patch(
                    master_geometry,
                    patch_id,
                    slave_position,
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

        const ID        best_patch_id = best_projection.patch_id;
        const Vec2      best_local    = best_projection.local;
        const Precision best_distance = best_projection.distance;
        const Precision gap           = best_projection.gap;

        diagnostics.maximum_closest_distance = std::max(
            diagnostics.maximum_closest_distance,
            best_distance
        );

        const NagataPatch& best_patch = master_geometry.patch(best_patch_id);
        const ProjectionInfo projection_info = classify_projection(best_projection);

        const auto slave_area_it = slave_tributary_areas.find(slave_node_id);
        const Precision slave_weight =
            slave_area_it != slave_tributary_areas.end()
                ? slave_area_it->second
                : Precision(1);

        const ContactLaw contact_law = evaluate_augmented_lagrange_law(
            gap,
            normal_multiplier,
            penalty
        );

        current_partners[slave_node_id] = best_patch_id;
        current_multipliers[slave_node_id] = normal_multiplier;
        current_gaps[slave_node_id] = gap;
        current_characteristic_lengths[slave_node_id] =
            std::sqrt(std::max(best_patch.area, Precision(0)));

        if (!contact_law.active) {
            ++diagnostics.open_closest_partner;
            continue;
        }

        ++diagnostics.active_contacts;
        current_active_slaves.insert(slave_node_id);

        const Precision penetration = std::max(Precision(0), -gap);
        const Precision slave_force = slave_weight * contact_law.pressure;

        diagnostics.slave_force_squared_sum += slave_force * slave_force;
        diagnostics.maximum_slave_force = std::max(
            diagnostics.maximum_slave_force,
            slave_force
        );

        if (penetration > diagnostics.maximum_penetration) {
            diagnostics.maximum_penetration = penetration;
            diagnostics.worst_slave         = slave_node_id;
            diagnostics.worst_surface       = best_patch_id;
            diagnostics.worst_local         = best_local;
            diagnostics.worst_mode          = projection_info.mode;
            diagnostics.worst_gap           = gap;
            diagnostics.worst_distance      = best_distance;
        }

        if (!previously_active) {
            ++diagnostics.activations;
        }

        if (previous_partner != state.partners.end() &&
            previous_partner->second != best_patch_id) {
            ++diagnostics.partner_switches;

            if (print_contact_details && detail_prints < maximum_detail_prints) {
                const auto previous_flags = std::cout.flags();
                const auto previous_precision = std::cout.precision();

                std::cout
                    << std::scientific
                    << std::setprecision(9)
                    << "[CONTACT_SWITCH]"
                    << " call="      << runtime_state.call
                    << " slave="     << slave_node_id
                    << " old="       << previous_partner->second
                    << " new="       << best_patch_id
                    << " new_mode="  << projection_mode_name(projection_info.mode)
                    << " new_dist="  << best_distance
                    << " new_gap="   << gap
                    << " new_local=(" << best_local(0)
                    << ','             << best_local(1)
                    << ")"
                    << '\n';

                std::cout.flags(previous_flags);
                std::cout.precision(previous_precision);
                ++detail_prints;
            }
        }

        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(slave_node_id))
        );

        hash_value(
            diagnostics.signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(best_patch_id))
        );

        assemble_contact_pair(
            master_geometry,
            best_patch,
            slave_node_id,
            best_local,
            projection_info,
            penalty,
            clearance,
            normal_multiplier,
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
        std::sqrt(diagnostics.slave_force_squared_sum);

    const auto assembly_end = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            assembly_end - assembly_start
        ).count();

    const std::size_t added_triplets = triplets.size() - initial_triplet_count;
    const Precision average_candidates =
        diagnostics.slave_nodes > 0
            ? static_cast<Precision>(diagnostics.candidate_surfaces) /
              static_cast<Precision>(diagnostics.slave_nodes)
            : Precision(0);

    if (print_contact_summary) {
        const auto previous_flags = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT]"
            << " call="          << runtime_state.call
            << " depth="         << runtime_state.trials.size()
            << " frozen="        << (freeze_surface_partners ? 1 : 0)
            << " patches="       << master_geometry.patch_count()
            << " active="        << diagnostics.active_contacts
            << " d_active="      << active_change
            << " activated="     << diagnostics.activations
            << " deactivated="   << diagnostics.deactivations
            << " switches="      << diagnostics.partner_switches
            << " changed="       << (signature_changed ? 1 : 0)
            << " lost_active="   << diagnostics.active_partner_losses
            << " retained_pen="  << diagnostics.retained_penetrations
            << " self_rej="      << diagnostics.self_contact_rejections
            << " interior="      << diagnostics.interior_consistent
            << " interior_fb="   << diagnostics.interior_fallback
            << " edge="          << diagnostics.edge_consistent
            << " edge_fb="       << diagnostics.edge_fallback
            << " corner="        << diagnostics.corner_direct
            << " K_full="        << diagnostics.consistent_tangents
            << " K_stab="        << diagnostics.stabilized_tangents
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
                << " patch="    << diagnostics.worst_surface
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
}

} // namespace fem::constraint
