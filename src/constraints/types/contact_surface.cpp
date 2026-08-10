/**
 * @file contact_surface.cpp
 * @brief Implements frictionless dual-mortar augmented-Lagrange contact.
 *
 * SurfaceRegion slaves use segment-to-segment integration over the projected
 * overlap of slave and master surface patches. Quadrature points are numerical
 * integration points only: they do not own contact state, active-set flags or
 * augmented-Lagrange multipliers.
 *
 * A local dual basis is constructed on every complete slave element from its
 * surface mass matrix. The normal gap is integrated against that dual basis and
 * accumulated into one normalized mortar constraint per global slave node.
 * Augmented multipliers and the unilateral active set therefore live exclusively
 * on those mortar constraints. Contact residuals and a symmetric frozen-geometry
 * tangent are assembled directly from the integrated mortar coupling vectors.
 *
 * Master/slave overlap is recomputed for every nonlinear evaluation. No master
 * facet is frozen for mortar contact; trial states are retained only for clean
 * commit/rollback semantics during predictors, line searches and cutbacks.
 * Explicit NodeRegion slaves continue to use Contact::assemble() in contact.cpp.
 *
 * @author Finn Eggers
 * @date 10.08.2026
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

constexpr Precision augmentation_gap_relative_tolerance        = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance        = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

constexpr bool print_contact_summary = true;

using LocalTriangle = std::array<Vec2, 3>;

struct MasterPatch {
    ID                  surface_id = -1;
    LocalTriangle       local {};
    std::array<Vec3, 3> global {};
    BvhAabb::Aabb       box = BvhAabb::Aabb::invalid();
};

struct MortarConstraintData {
    // Complete slave support used to normalize the dual weighted gap. This is
    // independent of the current overlap and therefore remains well conditioned
    // when contact starts on only a small fraction of one element.
    Precision support               = Precision(0);
    Precision overlap_measure       = Precision(0);
    Precision gap_integral          = Precision(0);
    Precision characteristic_length = Precision(0);

    // Derivative of the unnormalized weighted gap with frozen overlap geometry:
    //     H_i,a = int_Gamma phi_i B_a^T n dA.
    // The map key is one translational FE node participating in the constraint.
    std::unordered_map<ID, Vec3> gradient;
};

struct MortarDiagnostics {
    Index slave_surfaces     = 0;
    Index candidate_patches  = 0;
    Index maximum_candidates = 0;
    Index overlap_segments   = 0;
    Index quadrature_points  = 0;
    Index constraints        = 0;
    Index active_constraints = 0;
    Index activations        = 0;
    Index deactivations      = 0;
    Index self_rejections    = 0;
    Index invalid_dual_bases = 0;

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

    // Quadratic triangles are segmented at their midside nodes. The actual
    // curved FE mapping is still evaluated at every mortar quadrature point.
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

    // Serendipity quadrilaterals are segmented around the natural centre.
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
            if (std::abs(value) <= tangent_tolerance) {
                continue;
            }

            triplets.emplace_back(global_row, global_col, value);
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

/**
 * Builds a dimensionless dual multiplier basis on one complete slave element.
 *
 * With the primal surface mass matrix
 *
 *     M = int_Gamma N N^T dA
 *
 * and the positive diagonal scaling D_ii = M_ii, the transformation
 *
 *     A = D M^-1,   Phi = A N
 *
 * satisfies the local biorthogonality relation
 *
 *     int_Gamma Phi N^T dA = D.
 *
 * The matrix is formed on the complete slave element, not on the current
 * contact overlap. Partial contact can therefore never make this projection
 * singular. `support` returns diag(D) and is used to normalize the global nodal
 * mortar gap to units of length.
 */
bool build_dual_basis(
    const model::SurfaceInterface::Ptr&       surface,
    const model::Field&                       node_coords,
    const math::quadrature::Quadrature&       quadrature,
    DynamicMatrix&                            dual,
    DynamicVector&                            support
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

    support = mass.diagonal();
    dual    = inverse;

    for (Index i = 0; i < n; ++i) {
        if (!(support(i) > Precision(1e-14) * scale) || !std::isfinite(support(i))) {
            return false;
        }
        dual.row(i) *= support(i);
    }

    return dual.allFinite();
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

/**
 * Assembles frictionless dual-mortar contact for a SurfaceRegion slave.
 *
 * The algorithm has no persistent quadrature-point contact state:
 *
 *  1. build a dual basis on each complete slave element,
 *  2. integrate dual weighted gaps and coupling vectors over all projected
 *     slave/master overlap segments,
 *  3. normalize the accumulated weighted gap into one constraint per global
 *     slave node,
 *  4. evaluate the unilateral augmented pressure at those constraints,
 *  5. assemble residual and the fixed-overlap/fixed-normal consistent tangent
 *     from the integrated coupling vectors.
 *
 * Active-set changes are handled directly by this semismooth nodal law during
 * Newton. They are not master-partner changes and do not request an outer
 * active-set restart. The overlap geometry is recomputed on every assembly.
 */
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

    const auto        assembly_start       = std::chrono::steady_clock::now();
    const std::size_t initial_triplet_count = triplets.size();

    ++runtime_state.call;

    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    const Precision search_radius = distance + std::max(clearance, Precision(0));

    // ---------------------------------------------------------------------
    // Triangulated master broadphase in the current configuration
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
    std::unordered_map<ID, MortarConstraintData> constraints;

    std::vector<ID> candidate_patch_ids;
    candidate_patch_ids.reserve(64);

    // ---------------------------------------------------------------------
    // Build dual bases and integrate all projected overlap segments
    // ---------------------------------------------------------------------
    for (ID slave_surface_id : *slave_surfaces) {
        if (!valid_surface_id(slave_surface_id, surfaces)) {
            continue;
        }

        const auto& slave = surfaces[static_cast<std::size_t>(slave_surface_id)];
        ++diagnostics.slave_surfaces;

        DynamicMatrix dual;
        DynamicVector local_support;

        if (!build_dual_basis(
                slave,
                node_coords,
                triangle_quadrature,
                dual,
                local_support)) {
            ++diagnostics.invalid_dual_bases;
            logging::warning(false,
                "CONTACT: failed to construct dual mortar basis for slave surface ",
                slave_surface_id);
            continue;
        }

        const Precision slave_area   = std::max(slave->area(node_coords), Precision(0));
        const Precision slave_length = std::sqrt(slave_area);

        // The normalization support is assembled on the complete slave surface,
        // independently of how much of the element is currently in contact.
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

        for (ID patch_id : candidates) {
            if (patch_id < 0 || static_cast<std::size_t>(patch_id) >= master_patches.size()) {
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

            if (!overlap.is_ccw()) {
                overlap.flip();
            }

            const Index quadrature_before = diagnostics.quadrature_points;

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
                    const DynamicVector Phi = dual * Ns;

                    const Precision gap = (global - master_global).dot(normal) - clearance;

                    if (!Ns.allFinite() || !Nm.allFinite() || !Phi.allFinite() ||
                        !std::isfinite(gap) || !std::isfinite(weight) || weight <= Precision(0)) {
                        return;
                    }

                    ++diagnostics.quadrature_points;
                    diagnostics.maximum_geometric_penetration =
                        std::max(diagnostics.maximum_geometric_penetration,
                                 std::max(Precision(0), -gap));

                    // Each dual basis function contributes to one nodal mortar
                    // constraint. Quadrature points themselves own no state.
                    for (Index constraint_node = 0; constraint_node < slave->n_nodes; ++constraint_node) {
                        const Precision phi = Phi(constraint_node);
                        if (std::abs(phi) <= Precision(1e-14)) {
                            continue;
                        }

                        const ID constraint_id = slave->nodes()[constraint_node];
                        MortarConstraintData& constraint = constraints[constraint_id];
                        const Precision factor = weight * phi;

                        constraint.overlap_measure += weight * std::abs(phi);
                        constraint.gap_integral    += factor * gap;

                        for (Index local_node = 0; local_node < slave->n_nodes; ++local_node) {
                            add_gradient(
                                constraint.gradient,
                                slave->nodes()[local_node],
                                factor * Ns(local_node) * normal
                            );
                        }

                        for (Index local_node = 0; local_node < master->n_nodes; ++local_node) {
                            add_gradient(
                                constraint.gradient,
                                master->nodes()[local_node],
                                -factor * Nm(local_node) * normal
                            );
                        }
                    }
                }
            );

            if (diagnostics.quadrature_points > quadrature_before) {
                ++diagnostics.overlap_segments;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Evaluate nodal mortar constraints and assemble AL residual/tangent
    // ---------------------------------------------------------------------
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

        current_gaps[constraint_id] = gap;
        current_characteristic_lengths[constraint_id] =
            constraint.characteristic_length;

        const auto multiplier_it = state.normal_multipliers.find(constraint_id);
        const Precision normal_multiplier =
            multiplier_it != state.normal_multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        // One semismooth unilateral law per mortar constraint. No quadrature
        // point is activated or deactivated independently.
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

        // f_c = -p_i H_i. `gradient` already contains the complete mortar
        // integral over every overlap segment contributing to constraint i.
        for (const auto& [node_id, gradient] : constraint.gradient) {
            const Vec3 force = -pressure * gradient;
            add_translational_force(nodal_forces, node_id, force);
            diagnostics.force_contribution_squared_sum += force.squaredNorm();
        }

        // With frozen normal, overlap and dual basis during one linearization,
        //
        //   dg_i = H_i^T du / support_i,
        //   dp_i = -penalty dg_i,
        //   df_c  = penalty/support_i H_i H_i^T du.
        //
        // This is the consistent symmetric tangent of the discrete nodal mortar
        // law used for the residual above; it is not the old pointwise penalty
        // tangent.
        const Precision tangent_scale = penalty / constraint.support;

        for (const auto& [row_node, row_gradient] : constraint.gradient) {
            for (const auto& [col_node, col_gradient] : constraint.gradient) {
                const Mat3 block =
                    tangent_scale * row_gradient * col_gradient.transpose();

                add_tangent_block(
                    row_node,
                    col_node,
                    block,
                    system_nodal_dofs,
                    triplets
                );
            }
        }
    }

    // ---------------------------------------------------------------------
    // Update only geometric mortar state; AL multipliers stay fixed in Newton
    // ---------------------------------------------------------------------
    diagnostics.activations = 0;
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
        current_active_slaves.begin(),
        current_active_slaves.end()
    );
    std::sort(active_node_ids.begin(), active_node_ids.end());

    std::uint64_t signature = UINT64_C(1469598103934665603);
    for (ID node_id : active_node_ids) {
        hash_value(signature,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(node_id)));
    }

    // Surface mortar has no discrete master ownership. The nodal active set is
    // part of the semismooth Newton law and therefore must not request the outer
    // partner-update restart used by legacy node-to-surface contact.
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
    const Precision force_norm =
        std::sqrt(diagnostics.force_contribution_squared_sum);
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
            << " call="        << runtime_state.call
            << " mode=MORTAR"
            << " depth="       << runtime_state.trials.size()
            << " surfaces="    << diagnostics.slave_surfaces
            << " segments="    << diagnostics.overlap_segments
            << " qpts="        << diagnostics.quadrature_points
            << " constraints=" << diagnostics.constraints
            << " active="      << diagnostics.active_constraints
            << " activated="   << diagnostics.activations
            << " deactivated=" << diagnostics.deactivations
            << " active_changed=" << (active_changed ? 1 : 0)
            << " cand_avg="    << average_candidates
            << " cand_max="    << diagnostics.maximum_candidates
            << " self_rej="    << diagnostics.self_rejections
            << " dual_fail="   << diagnostics.invalid_dual_bases
            << " max_geom_pen="   << diagnostics.maximum_geometric_penetration
            << " max_mortar_pen=" << diagnostics.maximum_mortar_penetration
            << " max_pcoef="      << diagnostics.maximum_pressure_coefficient
            << " force_norm="  << force_norm
            << " triplets="    << added_triplets
            << " ms="          << elapsed_ms
            << '\n';

        std::cout.flags(previous_flags);
        std::cout.precision(previous_precision);
    }
}

/**
 * Updates augmented multipliers for nodal mortar constraints after a converged
 * inner Newton solve.
 *
 * Unlike node-to-surface contact, mortar has no discrete partner update to
 * defer. The multiplier update therefore acts directly on the current normalized
 * mortar gaps. Multipliers belonging to constraints that no longer overlap the
 * master interface are removed without requesting a restart because those
 * multipliers already contribute no residual in the converged geometry.
 */
bool Contact::update_augmented_lagrange_surface() const {
    AssemblyState& state =
        runtime_state.trials.empty()
            ? runtime_state.committed
            : runtime_state.trials.back().state;

    // Remove history for constraints that are no longer represented by any
    // current slave/master overlap. This does not change the assembled residual.
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
        maximum_change =
            std::max(maximum_change, multiplier_change);

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
