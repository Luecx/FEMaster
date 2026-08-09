/**
 * @file contact_surface.cpp
 * @brief Implements frictionless surface-to-surface mortar contact.
 *
 * SurfaceRegion slaves use segment-to-segment integration over the projected
 * overlap of slave and master surface patches. The geometric integration is
 * built on Surface::integrate_triangular(), which clips a projected master
 * triangle against the slave natural domain and returns the complete physical
 * integration weight.
 *
 * The scalar normal gap is projected to the slave nodes with the local mortar
 * mass matrix. Augmented-Lagrange multipliers are therefore stored per slave
 * node instead of per quadrature point. Contact residuals and a symmetric
 * frozen-normal tangent are integrated over the complete overlap. Explicit
 * NodeRegion slaves continue to use Contact::assemble() in contact.cpp.
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

namespace fem {
namespace constraint {
namespace {

constexpr Precision geometry_tolerance = Precision(1e-12);
constexpr Precision tangent_tolerance  = Precision(1e-14);
constexpr bool      print_contact_summary = true;

using LocalTriangle = std::array<Vec2, 3>;

struct MasterPatch {
    ID                  surface_id = -1;
    LocalTriangle       local {};
    std::array<Vec3, 3> global {};
    BvhAabb::Aabb       box = BvhAabb::Aabb::invalid();
};

struct MortarPoint {
    ID            master_surface_id = -1;
    Vec2          slave_local = Vec2::Zero();
    Vec2          master_local = Vec2::Zero();
    DynamicVector slave_shape;
    DynamicVector master_shape;
    Vec3          normal = Vec3::Zero();
    Precision     gap    = Precision(0);
    Precision     weight = Precision(0);
};

struct MortarDiagnostics {
    Index slave_surfaces        = 0;
    Index candidate_patches     = 0;
    Index maximum_candidates    = 0;
    Index overlap_segments      = 0;
    Index integration_points    = 0;
    Index active_points         = 0;
    Index active_nodes          = 0;
    Index activations           = 0;
    Index deactivations         = 0;
    Index self_rejections       = 0;
    Index regularized_gap_solves = 0;

    Precision maximum_penetration = Precision(0);
    Precision maximum_pressure    = Precision(0);
    Precision point_force_squared_sum = Precision(0);
};

bool valid_surface_id(
    ID                                                        surface_id,
    const std::vector<model::SurfaceInterface::Ptr>&          surfaces
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

Vec3 barycentric(const Vec2& point, const LocalTriangle& triangle) {
    const Vec2 edge_r = triangle[1] - triangle[0];
    const Vec2 edge_s = triangle[2] - triangle[0];
    const Vec2 offset = point - triangle[0];

    const Precision denominator = cross_2d(edge_r, edge_s);
    const Precision lambda_1    = cross_2d(offset, edge_s) / denominator;
    const Precision lambda_2    = cross_2d(edge_r, offset) / denominator;

    Vec3 lambda;
    lambda << Precision(1) - lambda_1 - lambda_2, lambda_1, lambda_2;
    return lambda;
}

std::vector<LocalTriangle> local_triangles(
    const model::SurfaceInterface::Ptr& surface
) {
    std::vector<LocalTriangle> triangles;

    // Quadratic triangles are split at their midside nodes so the master patch
    // used for overlap detection follows the curved surface more closely.
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

    // Serendipity quadrilaterals are split around the natural centre. The true
    // surface mapping is still evaluated at every integration point; these
    // triangles are used only to define the overlap segmentation.
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

BvhAabb::Aabb make_surface_aabb(
    const model::SurfaceInterface::Ptr& surface,
    const model::Field&                 node_coords
) {
    BvhAabb::Aabb box = BvhAabb::Aabb::invalid();

    for (Index local_node = 0; local_node < surface->n_nodes; ++local_node) {
        const ID node_id = surface->nodes()[local_node];
        box.expand_point(node_coords.row_vec3(static_cast<Index>(node_id)));
    }

    return box;
}

MasterPatch make_master_patch(
    ID                                  surface_id,
    const model::SurfaceInterface::Ptr& surface,
    const LocalTriangle&                local,
    const model::Field&                 node_coords
) {
    MasterPatch patch;
    patch.surface_id = surface_id;
    patch.local      = local;

    for (std::size_t i = 0; i < patch.local.size(); ++i) {
        patch.global[i] = surface->local_to_global(patch.local[i], node_coords);
        patch.box.expand_point(patch.global[i]);
    }

    return patch;
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
    ID              row_node,
    ID              col_node,
    const Mat3&     block,
    SystemDofIds&   system_nodal_dofs,
    TripletList&    triplets
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
            if (std::abs(value) <= tangent_tolerance) {
                continue;
            }

            triplets.emplace_back(global_row, global_col, value);
        }
    }
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
    logging::error(
        slave_surfaces != nullptr,
        "CONTACT: mortar surface assembly requires a SurfaceRegion slave"
    );

    logging::error(
        slave_nodes == nullptr,
        "CONTACT: mortar surface assembly cannot be used with a NodeRegion slave"
    );

    logging::error(
        model_data.positions != nullptr,
        "CONTACT: positions field not set in model data"
    );

    const model::Field& node_coords = *model_data.positions;
    auto&               surfaces    = model_data.surfaces;

    const auto assembly_start = std::chrono::steady_clock::now();
    const std::size_t initial_triplet_count = triplets.size();

    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    const bool frozen_trial =
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_surface_partners;

    const Precision search_radius =
        distance + std::max(clearance, Precision(0));

    // ---------------------------------------------------------------------
    // Build a triangulated master patch BVH in the current configuration
    // ---------------------------------------------------------------------
    std::vector<MasterPatch> master_patches;
    BvhAabb                  master_bvh;

    for (ID master_surface_id : *master_surfaces) {
        if (!valid_surface_id(master_surface_id, surfaces)) {
            continue;
        }

        const auto& master = surfaces[static_cast<std::size_t>(master_surface_id)];

        for (const auto& local_triangle : local_triangles(master)) {
            const ID patch_id = static_cast<ID>(master_patches.size());
            master_patches.push_back(
                make_master_patch(master_surface_id, master, local_triangle, node_coords)
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

    std::unordered_map<ID, Precision> gap_sum;
    std::unordered_map<ID, Precision> gap_weight;
    std::unordered_map<ID, Precision> characteristic_length;
    std::unordered_map<ID, Precision> representative_weight;
    std::unordered_map<ID, ID>        representative_master;

    std::vector<ID> candidate_patch_ids;
    candidate_patch_ids.reserve(64);

    // ---------------------------------------------------------------------
    // Integrate every slave/master overlap and assemble the mortar residual
    // ---------------------------------------------------------------------
    for (ID slave_surface_id : *slave_surfaces) {
        if (!valid_surface_id(slave_surface_id, surfaces)) {
            continue;
        }

        const auto& slave = surfaces[static_cast<std::size_t>(slave_surface_id)];
        ++diagnostics.slave_surfaces;

        BvhAabb::Aabb slave_box = make_surface_aabb(slave, node_coords);
        slave_box.inflate(search_radius);

        const auto& candidates = master_bvh.query_aabb(slave_box, &candidate_patch_ids);
        const Index candidate_count = static_cast<Index>(candidates.size());

        diagnostics.candidate_patches += candidate_count;
        diagnostics.maximum_candidates =
            std::max(diagnostics.maximum_candidates, candidate_count);

        const Index n_slave = slave->n_nodes;
        DynamicMatrix G = DynamicMatrix::Zero(n_slave, n_slave);
        DynamicVector h = DynamicVector::Zero(n_slave);

        std::vector<MortarPoint> mortar_points;
        mortar_points.reserve(static_cast<std::size_t>(std::max<Index>(candidate_count, 1)));

        for (ID patch_id : candidates) {
            if (patch_id < 0 ||
                static_cast<std::size_t>(patch_id) >= master_patches.size()) {
                continue;
            }

            const MasterPatch& patch = master_patches[static_cast<std::size_t>(patch_id)];
            if (!valid_surface_id(patch.surface_id, surfaces)) {
                continue;
            }

            const auto& master = surfaces[static_cast<std::size_t>(patch.surface_id)];

            if (surfaces_share_node(*slave, *master)) {
                ++diagnostics.self_rejections;
                continue;
            }

            LocalTriangle slave_local {};
            for (std::size_t i = 0; i < slave_local.size(); ++i) {
                slave_local[i] = slave->global_to_local(patch.global[i], node_coords, false);
            }

            if (std::abs(cross_2d(slave_local[1] - slave_local[0],
                                  slave_local[2] - slave_local[0])) <= geometry_tolerance) {
                continue;
            }

            model::SurfaceInterface::Polygon overlap {
                slave_local[0], slave_local[1], slave_local[2]
            };

            // Opposing contact faces naturally reverse the projected triangle
            // orientation. SurfacePolygon clipping requires a CCW subject.
            if (!overlap.is_ccw()) {
                overlap.flip();
            }

            const std::size_t points_before = mortar_points.size();

            slave->integrate_triangular(
                node_coords,
                overlap,
                triangle_quadrature,
                [&](const Vec2& local, const Vec3& global, Precision weight) {
                    const Vec3 lambda = barycentric(local, slave_local);
                    const Vec2 master_local =
                        lambda(0) * patch.local[0] +
                        lambda(1) * patch.local[1] +
                        lambda(2) * patch.local[2];

                    const Vec3 master_global = master->local_to_global(master_local, node_coords);
                    if ((master_global - global).norm() > search_radius + geometry_tolerance) {
                        return;
                    }

                    Vec3 normal = master->normal(node_coords, master_local);
                    if (!normal.allFinite() || normal.norm() <= geometry_tolerance) {
                        return;
                    }
                    normal.normalize();
                    if (flip_normal) {
                        normal = -normal;
                    }

                    const DynamicVector Ns = slave->shape_function(local);
                    const DynamicVector Nm = master->shape_function(master_local);
                    const Precision gap = (global - master_global).dot(normal) - clearance;

                    if (!Ns.allFinite() || !Nm.allFinite() || !std::isfinite(gap) ||
                        !std::isfinite(weight) || weight <= Precision(0)) {
                        return;
                    }

                    MortarPoint point;
                    point.master_surface_id = patch.surface_id;
                    point.slave_local       = local;
                    point.master_local      = master_local;
                    point.slave_shape       = Ns;
                    point.master_shape      = Nm;
                    point.normal            = normal;
                    point.gap               = gap;
                    point.weight            = weight;
                    mortar_points.push_back(std::move(point));

                    for (Index row = 0; row < n_slave; ++row) {
                        h(row) += weight * Ns(row) * gap;

                        for (Index col = 0; col < n_slave; ++col) {
                            G(row, col) += weight * Ns(row) * Ns(col);
                        }
                    }
                }
            );

            if (mortar_points.size() > points_before) {
                ++diagnostics.overlap_segments;
            }
        }

        diagnostics.integration_points +=
            static_cast<Index>(mortar_points.size());

        if (mortar_points.empty()) {
            continue;
        }

        // -----------------------------------------------------------------
        // Project the scalar gap onto the slave interpolation space
        // -----------------------------------------------------------------
        Precision matrix_scale = Precision(0);
        for (Index i = 0; i < n_slave; ++i) {
            matrix_scale = std::max(matrix_scale, std::abs(G(i, i)));
        }

        if (matrix_scale <= std::numeric_limits<Precision>::epsilon()) {
            continue;
        }

        DynamicMatrix regularized_G = G;
        regularized_G.diagonal().array() += Precision(1e-10) * matrix_scale;

        const Eigen::LDLT<DynamicMatrix> decomposition(regularized_G);
        DynamicVector local_gap = DynamicVector::Zero(n_slave);

        bool projected_gap_valid = decomposition.info() == Eigen::Success;
        if (projected_gap_valid) {
            local_gap = decomposition.solve(h);
            projected_gap_valid = local_gap.allFinite();
        }

        if (!projected_gap_valid) {
            ++diagnostics.regularized_gap_solves;

            for (Index local_node = 0; local_node < n_slave; ++local_node) {
                const Precision diagonal = G(local_node, local_node);
                if (std::abs(diagonal) > Precision(1e-14) * matrix_scale) {
                    local_gap(local_node) = h(local_node) / diagonal;
                }
            }
        }

        const Precision slave_area = std::max(slave->area(node_coords), Precision(0));
        const Precision slave_length = std::sqrt(slave_area);
        const Precision support_tolerance =
            Precision(1e-12) * std::max(slave_area, Precision(1e-12));

        for (Index local_node = 0; local_node < n_slave; ++local_node) {
            const Precision support = G(local_node, local_node);
            if (!(support > support_tolerance)) {
                continue;
            }

            const ID node_id = slave->nodes()[local_node];
            gap_sum[node_id]    += support * local_gap(local_node);
            gap_weight[node_id] += support;
            characteristic_length[node_id] =
                std::max(characteristic_length[node_id], slave_length);
        }

        // Pick one representative master id per constrained slave node. The
        // value is retained only because the common AL state iterates the
        // `partners` map; mortar topology itself is represented by the active
        // slave-node signature and is not frozen to one master facet.
        for (const MortarPoint& point : mortar_points) {
            for (Index local_node = 0; local_node < n_slave; ++local_node) {
                const Precision contribution =
                    point.weight * point.slave_shape(local_node) * point.slave_shape(local_node);

                const ID node_id = slave->nodes()[local_node];
                if (contribution > representative_weight[node_id]) {
                    representative_weight[node_id] = contribution;
                    representative_master[node_id] = point.master_surface_id;
                }
            }
        }

        // -----------------------------------------------------------------
        // Assemble AL traction and symmetric frozen-normal mortar tangent
        // -----------------------------------------------------------------
        DynamicVector nodal_multiplier = DynamicVector::Zero(n_slave);
        for (Index local_node = 0; local_node < n_slave; ++local_node) {
            const ID node_id = slave->nodes()[local_node];
            const auto multiplier_it = state.normal_multipliers.find(node_id);

            if (multiplier_it != state.normal_multipliers.end()) {
                nodal_multiplier(local_node) =
                    std::max(multiplier_it->second, Precision(0));
            }
        }

        for (const MortarPoint& point : mortar_points) {
            const auto& master = surfaces[static_cast<std::size_t>(point.master_surface_id)];

            Precision interpolated_multiplier =
                point.slave_shape.dot(nodal_multiplier);
            interpolated_multiplier =
                std::max(interpolated_multiplier, Precision(0));

            const Precision shifted_gap =
                point.gap - interpolated_multiplier / penalty;

            if (shifted_gap >= Precision(0)) {
                continue;
            }

            ++diagnostics.active_points;

            const Precision pressure = -penalty * shifted_gap;
            const Vec3 point_force =
                point.weight * penalty * shifted_gap * point.normal;

            diagnostics.maximum_penetration =
                std::max(diagnostics.maximum_penetration,
                         std::max(Precision(0), -point.gap));
            diagnostics.maximum_pressure =
                std::max(diagnostics.maximum_pressure, pressure);
            diagnostics.point_force_squared_sum +=
                point_force.squaredNorm();

            for (Index local_node = 0; local_node < n_slave; ++local_node) {
                add_translational_force(
                    nodal_forces,
                    slave->nodes()[local_node],
                    point.slave_shape(local_node) * point_force
                );
            }

            for (Index local_node = 0; local_node < master->n_nodes; ++local_node) {
                add_translational_force(
                    nodal_forces,
                    master->nodes()[local_node],
                    -point.master_shape(local_node) * point_force
                );
            }

            const Mat3 normal_stiffness =
                point.weight * penalty *
                point.normal * point.normal.transpose();

            for (Index row = 0; row < n_slave; ++row) {
                for (Index col = 0; col < n_slave; ++col) {
                    add_tangent_block(
                        slave->nodes()[row],
                        slave->nodes()[col],
                        point.slave_shape(row) * point.slave_shape(col) * normal_stiffness,
                        system_nodal_dofs,
                        triplets
                    );
                }

                for (Index col = 0; col < master->n_nodes; ++col) {
                    add_tangent_block(
                        slave->nodes()[row],
                        master->nodes()[col],
                        -point.slave_shape(row) * point.master_shape(col) * normal_stiffness,
                        system_nodal_dofs,
                        triplets
                    );
                }
            }

            for (Index row = 0; row < master->n_nodes; ++row) {
                for (Index col = 0; col < n_slave; ++col) {
                    add_tangent_block(
                        master->nodes()[row],
                        slave->nodes()[col],
                        -point.master_shape(row) * point.slave_shape(col) * normal_stiffness,
                        system_nodal_dofs,
                        triplets
                    );
                }

                for (Index col = 0; col < master->n_nodes; ++col) {
                    add_tangent_block(
                        master->nodes()[row],
                        master->nodes()[col],
                        point.master_shape(row) * point.master_shape(col) * normal_stiffness,
                        system_nodal_dofs,
                        triplets
                    );
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Convert element-local projected gaps into global nodal mortar state
    // ---------------------------------------------------------------------
    std::unordered_map<ID, ID>        current_partners;
    std::unordered_set<ID>            current_active_slaves;
    std::unordered_map<ID, Precision> current_multipliers;
    std::unordered_map<ID, Precision> current_gaps;
    std::unordered_map<ID, Precision> current_characteristic_lengths;

    current_partners.reserve(gap_sum.size());
    current_active_slaves.reserve(gap_sum.size());
    current_multipliers.reserve(gap_sum.size());
    current_gaps.reserve(gap_sum.size());
    current_characteristic_lengths.reserve(gap_sum.size());

    for (const auto& [node_id, weighted_gap] : gap_sum) {
        const Precision support = gap_weight[node_id];
        if (!(support > Precision(0))) {
            continue;
        }

        const auto master_it = representative_master.find(node_id);
        if (master_it == representative_master.end()) {
            continue;
        }

        const Precision gap = weighted_gap / support;
        const auto multiplier_it = state.normal_multipliers.find(node_id);
        const Precision normal_multiplier =
            multiplier_it != state.normal_multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        current_partners[node_id] = master_it->second;
        current_multipliers[node_id] = normal_multiplier;
        current_gaps[node_id] = gap;
        current_characteristic_lengths[node_id] = characteristic_length[node_id];

        if (gap - normal_multiplier / penalty < Precision(0)) {
            current_active_slaves.insert(node_id);
        }
    }

    diagnostics.active_nodes =
        static_cast<Index>(current_active_slaves.size());

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

    // Mortar integration may smoothly transfer area from one master facet to
    // its neighbour. Only the nodal unilateral active set is therefore treated
    // as a discrete nonlinear-state change.
    std::vector<ID> active_node_ids(
        current_active_slaves.begin(),
        current_active_slaves.end()
    );
    std::sort(active_node_ids.begin(), active_node_ids.end());

    std::uint64_t signature = UINT64_C(1469598103934665603);
    for (ID node_id : active_node_ids) {
        hash_value(
            signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(node_id))
        );
    }

    const bool signature_changed =
        state.previous_signature != 0 &&
        state.previous_signature != signature;

    state.partners                = std::move(current_partners);
    state.active_slaves           = std::move(current_active_slaves);
    state.normal_multipliers      = std::move(current_multipliers);
    state.gaps                    = std::move(current_gaps);
    state.characteristic_lengths = std::move(current_characteristic_lengths);
    state.previous_signature      = signature;
    state.previous_active         = diagnostics.active_nodes;
    state.last_signature_changed  = signature_changed;

    if (!frozen_trial &&
        !runtime_state.trials.empty() &&
        runtime_state.trials.back().freeze_after_update) {
        runtime_state.trials.back().freeze_surface_partners = true;
        runtime_state.trials.back().freeze_after_update     = false;
    }

    const auto assembly_end = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(assembly_end - assembly_start).count();

    const std::size_t added_triplets = triplets.size() - initial_triplet_count;
    const Precision force_norm = std::sqrt(diagnostics.point_force_squared_sum);
    const Precision average_candidates =
        diagnostics.slave_surfaces > 0
            ? static_cast<Precision>(diagnostics.candidate_patches) /
              static_cast<Precision>(diagnostics.slave_surfaces)
            : Precision(0);

    if (print_contact_summary) {
        const auto previous_flags = std::cout.flags();
        const auto previous_precision = std::cout.precision();

        std::cout
            << std::scientific
            << std::setprecision(3)
            << "[CONTACT]"
            << " call="       << runtime_state.call
            << " mode=MORTAR"
            << " depth="      << runtime_state.trials.size()
            << " frozen="     << (frozen_trial ? 1 : 0)
            << " surfaces="   << diagnostics.slave_surfaces
            << " segments="   << diagnostics.overlap_segments
            << " points="     << diagnostics.integration_points
            << " active="     << diagnostics.active_nodes
            << " active_qp="  << diagnostics.active_points
            << " activated="  << diagnostics.activations
            << " deactivated="<< diagnostics.deactivations
            << " cand_avg="   << average_candidates
            << " cand_max="   << diagnostics.maximum_candidates
            << " self_rej="   << diagnostics.self_rejections
            << " gap_reg="    << diagnostics.regularized_gap_solves
            << " max_pen="    << diagnostics.maximum_penetration
            << " max_press="  << diagnostics.maximum_pressure
            << " force_norm=" << force_norm
            << " changed="    << (signature_changed ? 1 : 0)
            << " triplets="   << added_triplets
            << " ms="         << elapsed_ms
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }
}

} // namespace constraint
} // namespace fem
