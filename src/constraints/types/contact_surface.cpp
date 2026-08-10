/**
 * @file contact_surface.cpp
 * @brief Implements frictionless dual-mortar augmented-Lagrange contact.
 *
 * SurfaceRegion slaves use a segment-to-segment mortar formulation. Both slave
 * and master subfacets are projected onto a common physical plane, clipped in
 * that plane and integrated over the resulting overlap polygons. Quadrature
 * points are numerical integration points only and own no contact history.
 *
 * A local dual basis is constructed on every complete slave element. Normal
 * gaps and mortar coupling vectors are integrated over the selected physical
 * overlap with exactly one visible master layer at every quadrature point.
 * Augmented multipliers and the unilateral active set live only on nodal mortar
 * constraints.
 *
 * DISTANCE is used only to bound the geometric search. It never activates
 * contact and it does not enter the signed normal gap. Surface mortar has no
 * persistent master-facet ownership and therefore no partner freezing.
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "contact.h"

#include "../../core/logging.h"
#include "../../model/geometry/surface/surface_polygon.h"
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

namespace fem {
namespace constraint {
namespace {

constexpr Precision geometry_tolerance = Precision(1e-12);
constexpr Precision area_tolerance     = Precision(1e-16);
constexpr Precision tangent_tolerance  = Precision(1e-14);

// Facing solid surfaces have opposite outward normals. A weak threshold rejects
// side/back faces while still allowing strongly curved opposing interfaces.
constexpr Precision maximum_opposing_normal_dot = Precision(-0.1);
constexpr Precision projection_selection_tolerance = Precision(1e-10);

constexpr Precision augmentation_gap_relative_tolerance        = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance        = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

constexpr bool print_contact_summary = true;

using LocalTriangle = std::array<Vec2, 3>;
using PlaneTriangle = std::array<Vec2, 3>;

struct SurfacePatch {
    ID                  surface_id = -1;
    LocalTriangle       local {};
    std::array<Vec3, 3> global {};
    BvhAabb::Aabb       box = BvhAabb::Aabb::invalid();
};

struct CommonPlane {
    Vec3 origin  = Vec3::Zero();
    Vec3 normal  = Vec3::Zero();
    Vec3 tangent = Vec3::Zero();
    Vec3 binormal = Vec3::Zero();
    bool valid = false;

    [[nodiscard]] Vec2 project(const Vec3& point) const {
        const Vec3 delta = point - origin;
        return Vec2(delta.dot(tangent), delta.dot(binormal));
    }
};

struct ProjectedMasterPatch {
    ID            patch_id = -1;
    PlaneTriangle plane {};
    SurfacePolygon<6> overlap;
};

struct SelectedMasterPoint {
    bool      valid      = false;
    ID        patch_id   = -1;
    ID        surface_id = -1;
    Vec2      local      = Vec2::Zero();
    Vec3      global     = Vec3::Zero();
    Vec3      normal     = Vec3::Zero();
    Precision distance   = std::numeric_limits<Precision>::max();
    Precision gap        = Precision(0);
};

struct MortarConstraintData {
    Precision support               = Precision(0);
    Precision overlap_measure       = Precision(0);
    Precision gap_integral          = Precision(0);
    Precision characteristic_length = Precision(0);
    Precision minimum_geometric_gap = std::numeric_limits<Precision>::max();

    std::unordered_map<ID, Vec3> gradient;
};

struct MortarDiagnostics {
    Index slave_surfaces       = 0;
    Index slave_patches        = 0;
    Index candidate_patches    = 0;
    Index maximum_candidates   = 0;
    Index projected_segments   = 0;
    Index overlap_segments     = 0;
    Index quadrature_points    = 0;
    Index constraints          = 0;
    Index active_constraints   = 0;
    Index activations          = 0;
    Index deactivations        = 0;
    Index self_rejections      = 0;
    Index normal_rejections    = 0;
    Index distance_rejections  = 0;
    Index hidden_layer_rejects = 0;
    Index invalid_dual_bases   = 0;

    Precision maximum_geometric_penetration = Precision(0);
    Precision maximum_mortar_penetration    = Precision(0);
    Precision maximum_pressure_coefficient  = Precision(0);
    Precision force_contribution_squared_sum = Precision(0);
};

bool valid_surface_id(
    ID                                               surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           static_cast<bool>(surfaces[static_cast<std::size_t>(surface_id)]);
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

Precision cross_2d(const Vec2& a, const Vec2& b) {
    return a(0) * b(1) - a(1) * b(0);
}

Vec3 barycentric(const Vec2& point, const PlaneTriangle& triangle) {
    const Vec2 edge_r = triangle[1] - triangle[0];
    const Vec2 edge_s = triangle[2] - triangle[0];
    const Vec2 offset = point - triangle[0];

    const Precision denominator = cross_2d(edge_r, edge_s);
    if (std::abs(denominator) <= geometry_tolerance) {
        return Vec3::Constant(std::numeric_limits<Precision>::quiet_NaN());
    }

    const Precision lambda_1 = cross_2d(offset, edge_s) / denominator;
    const Precision lambda_2 = cross_2d(edge_r, offset) / denominator;

    Vec3 lambda;
    lambda << Precision(1) - lambda_1 - lambda_2, lambda_1, lambda_2;
    return lambda;
}

bool barycentric_inside(const Vec3& lambda) {
    constexpr Precision tolerance = Precision(1e-10);
    return lambda.allFinite() &&
           lambda.minCoeff() >= -tolerance &&
           lambda.maxCoeff() <= Precision(1) + tolerance;
}

Vec2 interpolate_local(const LocalTriangle& triangle, const Vec3& lambda) {
    return lambda(0) * triangle[0] +
           lambda(1) * triangle[1] +
           lambda(2) * triangle[2];
}

std::vector<LocalTriangle> local_triangles(
    const model::SurfaceInterface::Ptr& surface
) {
    std::vector<LocalTriangle> triangles;

    if (surface->n_nodes == 6) {
        const Vec2 p0(Precision(0.0), Precision(0.0));
        const Vec2 p1(Precision(1.0), Precision(0.0));
        const Vec2 p2(Precision(0.0), Precision(1.0));
        const Vec2 p3(Precision(0.5), Precision(0.0));
        const Vec2 p4(Precision(0.5), Precision(0.5));
        const Vec2 p5(Precision(0.0), Precision(0.5));

        return {
            {p0, p3, p5},
            {p3, p1, p4},
            {p5, p4, p2},
            {p3, p4, p5}
        };
    }

    if (surface->n_nodes == 8) {
        const Vec2 center(Precision(0.0), Precision(0.0));
        const std::array<Vec2, 8> boundary {
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
}

SurfacePatch make_surface_patch(
    ID                                  surface_id,
    const model::SurfaceInterface::Ptr& surface,
    const LocalTriangle&                local,
    const model::Field&                 node_coords
) {
    SurfacePatch patch;
    patch.surface_id = surface_id;
    patch.local      = local;

    for (std::size_t i = 0; i < patch.local.size(); ++i) {
        patch.global[i] = surface->local_to_global(patch.local[i], node_coords);
        patch.box.expand_point(patch.global[i]);
    }
    return patch;
}

BvhAabb::Aabb make_surface_aabb(
    const model::SurfaceInterface::Ptr& surface,
    const model::Field&                 node_coords
) {
    BvhAabb::Aabb box = BvhAabb::Aabb::invalid();
    for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
        box.expand_point(node_coords.row_vec3(static_cast<Index>(surface->nodes()[local_node])));
    }
    return box;
}

CommonPlane make_common_plane(
    const model::SurfaceInterface::Ptr& slave,
    const SurfacePatch&                 slave_patch,
    const model::Field&                 node_coords
) {
    CommonPlane plane;

    const Vec2 local_center =
        (slave_patch.local[0] + slave_patch.local[1] + slave_patch.local[2]) /
        Precision(3);

    plane.origin = slave->local_to_global(local_center, node_coords);
    plane.normal = slave->normal(node_coords, local_center);

    if (!plane.origin.allFinite() || !plane.normal.allFinite() ||
        plane.normal.norm() <= geometry_tolerance) {
        return plane;
    }
    plane.normal.normalize();

    Vec3 first  = slave_patch.global[1] - slave_patch.global[0];
    Vec3 second = slave_patch.global[2] - slave_patch.global[0];
    first  -= first.dot(plane.normal)  * plane.normal;
    second -= second.dot(plane.normal) * plane.normal;

    plane.tangent = first.squaredNorm() >= second.squaredNorm() ? first : second;
    if (!plane.tangent.allFinite() || plane.tangent.norm() <= geometry_tolerance) {
        return plane;
    }
    plane.tangent.normalize();

    plane.binormal = plane.normal.cross(plane.tangent);
    if (!plane.binormal.allFinite() || plane.binormal.norm() <= geometry_tolerance) {
        return plane;
    }
    plane.binormal.normalize();
    plane.valid = true;
    return plane;
}

PlaneTriangle project_patch(const SurfacePatch& patch, const CommonPlane& plane) {
    return {
        plane.project(patch.global[0]),
        plane.project(patch.global[1]),
        plane.project(patch.global[2])
    };
}

Vec3 patch_normal(
    const SurfacePatch&                 patch,
    const model::SurfaceInterface::Ptr& surface,
    const model::Field&                 node_coords,
    bool                                flip
) {
    const Vec2 center =
        (patch.local[0] + patch.local[1] + patch.local[2]) / Precision(3);
    Vec3 normal = surface->normal(node_coords, center);
    if (flip) {
        normal = -normal;
    }
    if (normal.allFinite() && normal.norm() > geometry_tolerance) {
        normal.normalize();
    }
    return normal;
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

void add_tangent_block(
    ID            row_node,
    ID            col_node,
    const Mat3&   block,
    SystemDofIds& system_nodal_dofs,
    TripletList&  triplets
) {
    for (Dim local_row = 0; local_row < 3; ++local_row) {
        const int global_row = system_nodal_dofs(row_node, local_row);
        if (global_row < 0) {
            continue;
        }

        for (Dim local_col = 0; local_col < 3; ++local_col) {
            const int global_col = system_nodal_dofs(col_node, local_col);
            if (global_col < 0) {
                continue;
            }

            const Precision value = block(local_row, local_col);
            if (std::abs(value) > tangent_tolerance) {
                triplets.emplace_back(global_row, global_col, value);
            }
        }
    }
}

void add_gradient(
    std::unordered_map<ID, Vec3>& gradient,
    ID                            node_id,
    const Vec3&                   value
) {
    auto [it, inserted] = gradient.try_emplace(node_id, Vec3::Zero());
    (void) inserted;
    it->second += value;
}

bool build_dual_basis(
    const model::SurfaceInterface::Ptr& surface,
    const model::Field&                 node_coords,
    const math::quadrature::Quadrature& quadrature,
    DynamicMatrix&                      dual,
    DynamicVector&                      support
) {
    const Index n = surface->n_nodes;
    DynamicMatrix mass = DynamicMatrix::Zero(n, n);

    surface->integrate_triangular(
        node_coords,
        surface->local_domain_polygon(),
        quadrature,
        [&](const Vec2& local, const Vec3&, Precision weight) {
            const DynamicVector shape = surface->shape_function(local);
            mass.noalias() += weight * shape * shape.transpose();
        }
    );

    Precision scale = Precision(0);
    for (Index i = 0; i < n; ++i) {
        scale = std::max(scale, std::abs(mass(i, i)));
    }

    if (!(scale > std::numeric_limits<Precision>::epsilon()) || !mass.allFinite()) {
        return false;
    }

    const Eigen::LDLT<DynamicMatrix> decomposition(mass);
    if (decomposition.info() != Eigen::Success) {
        return false;
    }

    DynamicMatrix inverse = decomposition.solve(DynamicMatrix::Identity(n, n));
    if (!inverse.allFinite()) {
        return false;
    }

    // For linear S3/S4 surfaces the natural dual scaling is the row sum
    //
    //     D_i = sum_j M_ij = int_Gamma N_i dA,
    //
    // which gives sum_i D_i = area(Gamma). A constant multiplier therefore
    // reproduces a constant traction and passes the matching-surface patch test.
    // Quadratic S6/S8 shape functions can have zero or negative row integrals;
    // until a dedicated quadratic dual multiplier space is introduced, keep a
    // positive diagonal scaling but normalize it to the complete surface area.
    if (n == 3 || n == 4) {
        support = mass * DynamicVector::Ones(n);
    } else {
        support = mass.diagonal().cwiseAbs();

        const Precision surface_measure = mass.sum();
        const Precision support_measure = support.sum();
        if (!(surface_measure > Precision(0)) ||
            !(support_measure > Precision(0)) ||
            !std::isfinite(surface_measure) ||
            !std::isfinite(support_measure)) {
            return false;
        }
        support *= surface_measure / support_measure;
    }

    dual = inverse;

    for (Index i = 0; i < n; ++i) {
        if (!(support(i) > Precision(1e-14) * scale) || !std::isfinite(support(i))) {
            return false;
        }
        dual.row(i) *= support(i);
    }
    return dual.allFinite();
}

SelectedMasterPoint select_master_point(
    const Vec2&                                       plane_point,
    const Vec3&                                       slave_global,
    const Vec3&                                       slave_normal,
    const std::vector<ProjectedMasterPatch>&          projected_patches,
    const std::vector<SurfacePatch>&                  master_patches,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces,
    const model::Field&                               node_coords,
    Precision                                         search_radius,
    Precision                                         clearance,
    bool                                              flip_normal,
    MortarDiagnostics&                                diagnostics
) {
    SelectedMasterPoint best;

    for (const ProjectedMasterPatch& projected : projected_patches) {
        const Vec3 lambda = barycentric(plane_point, projected.plane);
        if (!barycentric_inside(lambda)) {
            continue;
        }

        if (projected.patch_id < 0 ||
            static_cast<std::size_t>(projected.patch_id) >= master_patches.size()) {
            continue;
        }

        const SurfacePatch& patch = master_patches[static_cast<std::size_t>(projected.patch_id)];
        if (!valid_surface_id(patch.surface_id, surfaces)) {
            continue;
        }

        const auto& master = surfaces[static_cast<std::size_t>(patch.surface_id)];
        const Vec2 master_local = interpolate_local(patch.local, lambda);
        const Vec3 master_global = master->local_to_global(master_local, node_coords);

        Vec3 master_normal = master->normal(node_coords, master_local);
        if (flip_normal) {
            master_normal = -master_normal;
        }

        if (!master_global.allFinite() || !master_normal.allFinite() ||
            master_normal.norm() <= geometry_tolerance) {
            continue;
        }
        master_normal.normalize();

        if (slave_normal.dot(master_normal) > maximum_opposing_normal_dot) {
            continue;
        }

        const Vec3 difference = slave_global - master_global;
        const Precision candidate_distance = difference.norm();
        if (!std::isfinite(candidate_distance) ||
            candidate_distance > search_radius + geometry_tolerance) {
            ++diagnostics.distance_rejections;
            continue;
        }

        const Precision candidate_gap = difference.dot(master_normal) - clearance;
        if (!std::isfinite(candidate_gap)) {
            continue;
        }

        const bool better =
            !best.valid ||
            candidate_distance < best.distance - projection_selection_tolerance ||
            (std::abs(candidate_distance - best.distance) <= projection_selection_tolerance &&
             projected.patch_id < best.patch_id);

        if (!better) {
            continue;
        }

        best.valid      = true;
        best.patch_id   = projected.patch_id;
        best.surface_id = patch.surface_id;
        best.local      = master_local;
        best.global     = master_global;
        best.normal     = master_normal;
        best.distance   = candidate_distance;
        best.gap        = candidate_gap;
    }

    return best;
}

void hash_value(std::uint64_t& hash, std::uint64_t value) {
    constexpr std::uint64_t prime = UINT64_C(1099511628211);
    for (int byte = 0; byte < 8; ++byte) {
        hash ^= value & UINT64_C(0xff);
        hash *= prime;
        value >>= 8;
    }
}

} // namespace

void Contact::assemble_surface(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList&      triplets
) const {
    logging::error(slave_surfaces != nullptr,
        "CONTACT: mortar surface assembly requires a SurfaceRegion slave");
    logging::error(slave_nodes == nullptr,
        "CONTACT: mortar surface assembly cannot be used with a NodeRegion slave");
    logging::error(model_data.positions != nullptr,
        "CONTACT: positions field not set in model data");

    const model::Field& node_coords = *model_data.positions;
    auto&               surfaces    = model_data.surfaces;

    const auto assembly_start = std::chrono::steady_clock::now();
    const std::size_t initial_triplet_count = triplets.size();
    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    const Precision search_radius = distance + std::max(clearance, Precision(0));

    std::vector<SurfacePatch> master_patches;
    BvhAabb                  master_bvh;

    for (ID master_surface_id : *master_surfaces) {
        if (!valid_surface_id(master_surface_id, surfaces)) {
            continue;
        }

        const auto& master = surfaces[static_cast<std::size_t>(master_surface_id)];
        for (const auto& local_triangle : local_triangles(master)) {
            const ID patch_id = static_cast<ID>(master_patches.size());
            master_patches.push_back(
                make_surface_patch(master_surface_id, master, local_triangle, node_coords)
            );
            master_bvh.add_aabb(patch_id, master_patches.back().box);
        }
    }
    master_bvh.finalize();

    const math::quadrature::Quadrature triangle_quadrature {
        math::quadrature::DOMAIN_ISO_TRI,
        math::quadrature::ORDER_QUARTIC
    };

    MortarDiagnostics diagnostics;
    std::unordered_map<ID, MortarConstraintData> constraints;

    std::vector<ID> candidate_patch_ids;
    candidate_patch_ids.reserve(64);

    for (ID slave_surface_id : *slave_surfaces) {
        if (!valid_surface_id(slave_surface_id, surfaces)) {
            continue;
        }

        const auto& slave = surfaces[static_cast<std::size_t>(slave_surface_id)];
        ++diagnostics.slave_surfaces;

        DynamicMatrix dual;
        DynamicVector local_support;
        if (!build_dual_basis(slave, node_coords, triangle_quadrature, dual, local_support)) {
            ++diagnostics.invalid_dual_bases;
            logging::warning(false,
                "CONTACT: failed to construct dual mortar basis for slave surface ",
                slave_surface_id);
            continue;
        }

        const Precision slave_area   = std::max(slave->area(node_coords), Precision(0));
        const Precision slave_length = std::sqrt(slave_area);

        for (Index local_node = 0; local_node < slave->n_nodes; ++local_node) {
            const ID node_id = slave->nodes()[local_node];
            MortarConstraintData& constraint = constraints[node_id];
            constraint.support += local_support(local_node);
            constraint.characteristic_length =
                std::max(constraint.characteristic_length, slave_length);
        }

        if (!master_bvh.valid()) {
            continue;
        }

        BvhAabb::Aabb slave_box = make_surface_aabb(slave, node_coords);
        slave_box.inflate(search_radius);

        const auto& candidates = master_bvh.query_aabb(slave_box, &candidate_patch_ids);
        const Index candidate_count = static_cast<Index>(candidates.size());
        diagnostics.candidate_patches += candidate_count;
        diagnostics.maximum_candidates =
            std::max(diagnostics.maximum_candidates, candidate_count);

        for (const LocalTriangle& slave_local_triangle : local_triangles(slave)) {
            const SurfacePatch slave_patch =
                make_surface_patch(slave_surface_id, slave, slave_local_triangle, node_coords);
            ++diagnostics.slave_patches;

            const CommonPlane plane = make_common_plane(slave, slave_patch, node_coords);
            if (!plane.valid) {
                continue;
            }

            const PlaneTriangle slave_plane = project_patch(slave_patch, plane);
            if (std::abs(cross_2d(slave_plane[1] - slave_plane[0],
                                  slave_plane[2] - slave_plane[0])) <= area_tolerance) {
                continue;
            }

            SurfacePolygon<3> slave_polygon {
                slave_plane[0], slave_plane[1], slave_plane[2]
            };
            if (!slave_polygon.is_ccw()) {
                slave_polygon.flip();
            }

            std::vector<ProjectedMasterPatch> projected_patches;
            projected_patches.reserve(static_cast<std::size_t>(candidate_count));

            for (ID patch_id : candidates) {
                if (patch_id < 0 ||
                    static_cast<std::size_t>(patch_id) >= master_patches.size()) {
                    continue;
                }

                const SurfacePatch& master_patch =
                    master_patches[static_cast<std::size_t>(patch_id)];
                if (!valid_surface_id(master_patch.surface_id, surfaces)) {
                    continue;
                }

                const auto& master = surfaces[static_cast<std::size_t>(master_patch.surface_id)];
                if (surfaces_share_node(*slave, *master)) {
                    ++diagnostics.self_rejections;
                    continue;
                }

                const Vec3 master_center_normal =
                    patch_normal(master_patch, master, node_coords, flip_normal);
                if (!master_center_normal.allFinite() ||
                    master_center_normal.norm() <= geometry_tolerance ||
                    plane.normal.dot(master_center_normal) > maximum_opposing_normal_dot) {
                    ++diagnostics.normal_rejections;
                    continue;
                }

                const PlaneTriangle master_plane = project_patch(master_patch, plane);
                if (std::abs(cross_2d(master_plane[1] - master_plane[0],
                                      master_plane[2] - master_plane[0])) <= area_tolerance) {
                    continue;
                }

                SurfacePolygon<3> master_polygon {
                    master_plane[0], master_plane[1], master_plane[2]
                };
                if (!master_polygon.is_ccw()) {
                    master_polygon.flip();
                }

                const auto overlap = master_polygon.intersection(slave_polygon);
                if (overlap.size() < 3 || overlap.area() <= area_tolerance) {
                    continue;
                }

                ProjectedMasterPatch projected;
                projected.patch_id = patch_id;
                projected.plane    = master_plane;
                projected.overlap  = overlap;
                projected_patches.push_back(std::move(projected));
                ++diagnostics.projected_segments;
            }

            if (projected_patches.empty()) {
                continue;
            }

            for (const ProjectedMasterPatch& segment : projected_patches) {
                const Index selected_before = diagnostics.quadrature_points;
                const Vec2 origin = segment.overlap[0];

                for (std::size_t fan = 1; fan + 1 < segment.overlap.size(); ++fan) {
                    const Vec2 edge_r = segment.overlap[fan]     - origin;
                    const Vec2 edge_s = segment.overlap[fan + 1] - origin;
                    const Precision triangle_jacobian =
                        std::abs(cross_2d(edge_r, edge_s));

                    if (triangle_jacobian <= area_tolerance) {
                        continue;
                    }

                    for (Index q = 0; q < triangle_quadrature.count(); ++q) {
                        const auto quadrature_point = triangle_quadrature.get_point(q);
                        const Vec2 plane_point =
                            origin + quadrature_point.r * edge_r + quadrature_point.s * edge_s;
                        const Precision weight = quadrature_point.w * triangle_jacobian;

                        const Vec3 slave_lambda = barycentric(plane_point, slave_plane);
                        if (!barycentric_inside(slave_lambda)) {
                            continue;
                        }

                        const Vec2 slave_local =
                            interpolate_local(slave_patch.local, slave_lambda);
                        const Vec3 slave_global =
                            slave->local_to_global(slave_local, node_coords);
                        Vec3 slave_normal = slave->normal(node_coords, slave_local);

                        if (!slave_global.allFinite() || !slave_normal.allFinite() ||
                            slave_normal.norm() <= geometry_tolerance) {
                            continue;
                        }
                        slave_normal.normalize();

                        const SelectedMasterPoint selected = select_master_point(
                            plane_point,
                            slave_global,
                            slave_normal,
                            projected_patches,
                            master_patches,
                            surfaces,
                            node_coords,
                            search_radius,
                            clearance,
                            flip_normal,
                            diagnostics
                        );

                        if (!selected.valid || selected.patch_id != segment.patch_id) {
                            if (selected.valid) {
                                ++diagnostics.hidden_layer_rejects;
                            }
                            continue;
                        }

                        const auto& master =
                            surfaces[static_cast<std::size_t>(selected.surface_id)];
                        const DynamicVector Ns   = slave->shape_function(slave_local);
                        const DynamicVector Nm   = master->shape_function(selected.local);
                        const DynamicVector Phi  = dual * Ns;

                        if (!Ns.allFinite() || !Nm.allFinite() || !Phi.allFinite() ||
                            !std::isfinite(selected.gap) || !std::isfinite(weight) ||
                            weight <= Precision(0)) {
                            continue;
                        }

                        ++diagnostics.quadrature_points;
                        diagnostics.maximum_geometric_penetration =
                            std::max(diagnostics.maximum_geometric_penetration,
                                     std::max(Precision(0), -selected.gap));

                        for (Index constraint_node = 0;
                             constraint_node < slave->n_nodes;
                             ++constraint_node) {
                            const Precision phi = Phi(constraint_node);
                            if (std::abs(phi) <= Precision(1e-14)) {
                                continue;
                            }

                            const ID constraint_id = slave->nodes()[constraint_node];
                            MortarConstraintData& constraint = constraints[constraint_id];
                            const Precision factor = weight * phi;

                            constraint.overlap_measure += weight * std::abs(phi);
                            constraint.gap_integral    += factor * selected.gap;
                            constraint.minimum_geometric_gap =
                                std::min(constraint.minimum_geometric_gap, selected.gap);

                            for (Index local_node = 0; local_node < slave->n_nodes; ++local_node) {
                                add_gradient(
                                    constraint.gradient,
                                    slave->nodes()[local_node],
                                    factor * Ns(local_node) * selected.normal
                                );
                            }

                            for (Index local_node = 0; local_node < master->n_nodes; ++local_node) {
                                add_gradient(
                                    constraint.gradient,
                                    master->nodes()[local_node],
                                    -factor * Nm(local_node) * selected.normal
                                );
                            }
                        }
                    }
                }

                if (diagnostics.quadrature_points > selected_before) {
                    ++diagnostics.overlap_segments;
                }
            }
        }
    }

    std::unordered_set<ID>            current_active_slaves;
    std::unordered_map<ID, Precision> current_gaps;
    std::unordered_map<ID, Precision> current_characteristic_lengths;

    current_active_slaves.reserve(constraints.size());
    current_gaps.reserve(constraints.size());
    current_characteristic_lengths.reserve(constraints.size());

    for (auto& [constraint_id, constraint] : constraints) {
        const Precision support_tolerance =
            Precision(1e-14) * std::max(constraint.support, Precision(1));
        const Precision overlap_tolerance =
            Precision(1e-12) * std::max(constraint.support, Precision(1e-16));

        if (!(constraint.support > support_tolerance) ||
            !(constraint.overlap_measure > overlap_tolerance) ||
            constraint.gradient.empty()) {
            continue;
        }

        ++diagnostics.constraints;

        const Precision gap = constraint.gap_integral / constraint.support;
        if (!std::isfinite(gap)) {
            continue;
        }

        const auto multiplier_it = state.normal_multipliers.find(constraint_id);
        const Precision normal_multiplier =
            multiplier_it != state.normal_multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        const bool geometrically_closed =
            constraint.minimum_geometric_gap < -geometry_tolerance;
        const bool has_multiplier_history = normal_multiplier > Precision(0);

        if (!geometrically_closed && !has_multiplier_history) {
            continue;
        }

        current_gaps[constraint_id] = gap;
        current_characteristic_lengths[constraint_id] =
            constraint.characteristic_length;

        const Precision pressure =
            std::max(Precision(0), normal_multiplier - penalty * gap);

        diagnostics.maximum_mortar_penetration =
            std::max(diagnostics.maximum_mortar_penetration,
                     std::max(Precision(0), -gap));
        diagnostics.maximum_pressure_coefficient =
            std::max(diagnostics.maximum_pressure_coefficient, pressure);

        if (!(pressure > Precision(0))) {
            continue;
        }

        current_active_slaves.insert(constraint_id);
        ++diagnostics.active_constraints;

        for (const auto& [node_id, gradient] : constraint.gradient) {
            const Vec3 force = -pressure * gradient;
            add_translational_force(nodal_forces, node_id, force);
            diagnostics.force_contribution_squared_sum += force.squaredNorm();
        }

        const Precision tangent_scale = penalty / constraint.support;
        for (const auto& [row_node, row_gradient] : constraint.gradient) {
            for (const auto& [col_node, col_gradient] : constraint.gradient) {
                add_tangent_block(
                    row_node,
                    col_node,
                    tangent_scale * row_gradient * col_gradient.transpose(),
                    system_nodal_dofs,
                    triplets
                );
            }
        }
    }

    diagnostics.activations   = 0;
    diagnostics.deactivations = 0;

    for (ID node_id : current_active_slaves) {
        if (state.active_slaves.find(node_id) == state.active_slaves.end()) {
            ++diagnostics.activations;
        }
    }
    for (ID node_id : state.active_slaves) {
        if (current_active_slaves.find(node_id) == current_active_slaves.end()) {
            ++diagnostics.deactivations;
        }
    }

    const bool active_changed =
        diagnostics.activations > 0 || diagnostics.deactivations > 0;

    std::vector<ID> active_node_ids(
        current_active_slaves.begin(), current_active_slaves.end()
    );
    std::sort(active_node_ids.begin(), active_node_ids.end());

    std::uint64_t signature = UINT64_C(1469598103934665603);
    for (ID node_id : active_node_ids) {
        hash_value(signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(node_id)));
    }

    state.partners.clear();
    state.active_slaves           = std::move(current_active_slaves);
    state.gaps                    = std::move(current_gaps);
    state.characteristic_lengths = std::move(current_characteristic_lengths);
    state.previous_signature      = signature;
    state.previous_active         = diagnostics.active_constraints;
    state.last_signature_changed  = false;

    const auto assembly_end = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(assembly_end - assembly_start).count();

    const std::size_t added_triplets = triplets.size() - initial_triplet_count;
    const Precision force_norm = std::sqrt(diagnostics.force_contribution_squared_sum);
    const Precision average_candidates =
        diagnostics.slave_surfaces > 0
            ? static_cast<Precision>(diagnostics.candidate_patches) /
              static_cast<Precision>(diagnostics.slave_surfaces)
            : Precision(0);

    if (print_contact_summary) {
        const auto previous_flags     = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT]"
            << " call="          << runtime_state.call
            << " mode=MORTAR"
            << " depth="         << runtime_state.trials.size()
            << " surfaces="      << diagnostics.slave_surfaces
            << " spatches="      << diagnostics.slave_patches
            << " seg_raw="       << diagnostics.projected_segments
            << " segments="      << diagnostics.overlap_segments
            << " qpts="          << diagnostics.quadrature_points
            << " constraints="   << diagnostics.constraints
            << " active="        << diagnostics.active_constraints
            << " activated="     << diagnostics.activations
            << " deactivated="   << diagnostics.deactivations
            << " active_changed="<< (active_changed ? 1 : 0)
            << " cand_avg="      << average_candidates
            << " cand_max="      << diagnostics.maximum_candidates
            << " self_rej="      << diagnostics.self_rejections
            << " normal_rej="    << diagnostics.normal_rejections
            << " layer_rej="     << diagnostics.hidden_layer_rejects
            << " dist_rej="      << diagnostics.distance_rejections
            << " dual_fail="     << diagnostics.invalid_dual_bases
            << " max_geom_pen="  << diagnostics.maximum_geometric_penetration
            << " max_mortar_pen="<< diagnostics.maximum_mortar_penetration
            << " max_pcoef="     << diagnostics.maximum_pressure_coefficient
            << " force_norm="    << force_norm
            << " triplets="      << added_triplets
            << " ms="            << elapsed_ms
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }
}

bool Contact::update_augmented_lagrange_surface() const {
    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    for (auto it = state.normal_multipliers.begin();
         it != state.normal_multipliers.end();) {
        if (state.gaps.find(it->first) == state.gaps.end()) {
            it = state.normal_multipliers.erase(it);
        } else {
            ++it;
        }
    }

    Index     updated_constraints = 0;
    Precision maximum_penetration = Precision(0);
    Precision maximum_change      = Precision(0);

    for (const auto& [constraint_id, gap] : state.gaps) {
        if (!std::isfinite(gap)) {
            continue;
        }

        const auto length_it = state.characteristic_lengths.find(constraint_id);
        const Precision characteristic_length =
            length_it != state.characteristic_lengths.end()
                ? std::max(length_it->second, Precision(0))
                : Precision(0);

        const Precision gap_tolerance =
            std::max(
                augmentation_gap_absolute_tolerance,
                augmentation_gap_relative_tolerance * characteristic_length
            );

        const auto multiplier_it = state.normal_multipliers.find(constraint_id);
        const Precision old_multiplier =
            multiplier_it != state.normal_multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        Precision new_multiplier = old_multiplier;
        if (gap < -gap_tolerance ||
            (gap > gap_tolerance && old_multiplier > Precision(0))) {
            new_multiplier =
                std::max(Precision(0), old_multiplier - penalty * gap);
        }

        const Precision multiplier_change =
            std::abs(new_multiplier - old_multiplier);
        const Precision multiplier_scale =
            std::max({
                old_multiplier,
                new_multiplier,
                penalty * gap_tolerance,
                Precision(1)
            });
        const Precision multiplier_tolerance =
            std::max(
                augmentation_multiplier_absolute_tolerance,
                augmentation_multiplier_relative_tolerance * multiplier_scale
            );

        if (new_multiplier > Precision(0) || old_multiplier > Precision(0)) {
            state.normal_multipliers[constraint_id] = new_multiplier;
        }

        maximum_penetration =
            std::max(maximum_penetration, std::max(Precision(0), -gap));
        maximum_change = std::max(maximum_change, multiplier_change);

        if (multiplier_change > multiplier_tolerance) {
            ++updated_constraints;
        }
    }

    state.last_augmentation_changed = updated_constraints > 0;

    if (print_contact_summary) {
        const auto previous_flags     = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT_AUGMENT]"
            << " mode=MORTAR"
            << " updated="     << updated_constraints
            << " max_dlambda=" << maximum_change
            << " max_pen="     << maximum_penetration
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }

    return state.last_augmentation_changed;
}

} // namespace constraint
} // namespace fem