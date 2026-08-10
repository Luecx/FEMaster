/**
 * @file contact.cpp
 * @brief Implements frictionless dual-mortar surface-to-surface contact.
 *
 * The complete contact formulation is implemented in this translation unit.
 * Every current slave subfacet defines one tangent plane. The slave triangle and
 * every admissible master triangle are projected onto this common physical plane,
 * clipped there and integrated over the resulting overlap polygons. All master
 * surfaces are tested directly; no BVH, search radius or persistent master
 * ownership participates in the formulation.
 *
 * Slave surfaces are processed independently in OpenMP workers when OpenMP is
 * available. Every worker owns a private mortar-constraint map; the maps are
 * reduced after the geometric integration so no locks are required in the
 * quadrature-point hot path.
 *
 * A local dual basis is constructed on every complete slave element. The overlap
 * integration accumulates one normalized normal-gap constraint and its frozen-
 * geometry gradient on every global slave mortar node. The unilateral law is
 *
 *     p_i = max(0, lambda_i - epsilon g_i),
 *
 * and the corresponding residual and tangent contributions are
 *
 *     f_c = -sum_i p_i H_i,
 *     K_c =  sum_i epsilon / S_i H_i H_i^T.
 *
 * Augmented-Lagrange multipliers are the only persistent physical contact
 * history. Transactional copies support nonlinear increment and line-search
 * rollback while the geometric overlap is always recomputed.
 *
 * @see Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#include "contact.h"

#include "../../core/config.h"
#include "../../core/logging.h"
#include "../../math/quadrature.h"
#include "../../model/geometry/surface/surface_polygon.h"
#include "../../model/model_data.h"

#include <Eigen/Cholesky>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem::constraint {
namespace {

constexpr Precision geometry_tolerance = Precision(1e-12);
constexpr Precision area_tolerance     = Precision(1e-16);
constexpr Precision tangent_tolerance  = Precision(1e-14);

// Facing solid surfaces have opposite outward normals. The weak threshold keeps
// moderately curved opposing facets while rejecting side and back faces.
constexpr Precision maximum_opposing_normal_dot    = Precision(-0.1);
constexpr Precision projection_selection_tolerance = Precision(1e-10);

constexpr Precision augmentation_gap_relative_tolerance        = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance        = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

using LocalTriangle = std::array<Vec2, 3>;
using PlaneTriangle = std::array<Vec2, 3>;

/**
 * @brief One master subtriangle projected onto the active slave tangent plane.
 */
struct ProjectedSegment {
    ID                master_surface_id = ID(-1);
    LocalTriangle     master_local {};
    PlaneTriangle     master_plane {};
    SurfacePolygon<6> overlap;
};

/**
 * @brief Integrated data of one global slave mortar constraint.
 */
struct MortarConstraint {
    Precision support               = Precision(0);
    Precision overlap_measure       = Precision(0);
    Precision gap_integral          = Precision(0);
    Precision characteristic_length = Precision(0);
    Precision minimum_geometric_gap = std::numeric_limits<Precision>::max();

    std::unordered_map<ID, Vec3> gradient;
};

bool valid_surface_id(
    ID                                                surface_id,
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

        std::vector<LocalTriangle> triangles;
        triangles.reserve(boundary.size());

        for (std::size_t i = 0; i < boundary.size(); ++i) {
            triangles.push_back({
                boundary[i],
                boundary[(i + 1) % boundary.size()],
                center
            });
        }

        return triangles;
    }

    const auto domain = surface->local_domain_polygon();
    std::vector<LocalTriangle> triangles;

    if (domain.size() < 3) {
        return triangles;
    }

    triangles.reserve(domain.size() - 2);
    for (std::size_t i = 1; i + 1 < domain.size(); ++i) {
        triangles.push_back({domain[0], domain[i], domain[i + 1]});
    }

    return triangles;
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

    // Linear S3/S4 surfaces use D_i = int N_i dA. For quadratic S6/S8
    // surfaces keep the existing positive full-area normalization.
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

void add_gradient(
    std::unordered_map<ID, Vec3>& gradient,
    ID                            node_id,
    const Vec3&                   value
) {
    auto [it, inserted] = gradient.try_emplace(node_id, Vec3::Zero());
    (void) inserted;
    it->second += value;
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

Contact::State& Contact::state() const {
    return trial_states.empty() ? committed_state : trial_states.back();
}

void Contact::begin_trial() const {
    trial_states.push_back(state());
}

void Contact::commit_trial() const {
    logging::error(!trial_states.empty(), "CONTACT: no trial state to commit");

    State accepted = std::move(trial_states.back());
    trial_states.pop_back();

    if (trial_states.empty()) {
        committed_state = std::move(accepted);
    } else {
        trial_states.back() = std::move(accepted);
    }
}

void Contact::rollback_trial() const {
    logging::error(!trial_states.empty(), "CONTACT: no trial state to roll back");
    trial_states.pop_back();
}

void Contact::assemble(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList&      triplets
) const {
    // -------------------------------------------------------------------------
    // 1. Validate model state and initialize thread-local mortar accumulation
    // -------------------------------------------------------------------------
    logging::error(master_surfaces != nullptr && slave_surfaces != nullptr,
        "CONTACT: surface regions are not defined");
    logging::error(model_data.positions != nullptr,
        "CONTACT: POSITION field is not set");

    const model::Field& node_coords = *model_data.positions;
    auto&               surfaces    = model_data.surfaces;
    State&              current     = state();

    current.gaps.clear();
    current.characteristic_lengths.clear();

    const math::quadrature::Quadrature quadrature {
        math::quadrature::DOMAIN_ISO_TRI,
        math::quadrature::ORDER_QUARTIC
    };

    int num_threads = 1;
#ifdef _OPENMP
    num_threads = std::max(
        1,
        std::min(
            static_cast<int>(global_config.max_threads),
            static_cast<int>(slave_surfaces->size())
        )
    );
#endif

    std::vector<std::unordered_map<ID, MortarConstraint>> thread_constraints(
        static_cast<std::size_t>(num_threads)
    );

    std::atomic_bool   assembly_failed{false};
    std::exception_ptr assembly_exception = nullptr;

    // -------------------------------------------------------------------------
    // 2. Process slave surfaces independently
    // -------------------------------------------------------------------------
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
    for (Index slave_index = 0;
         slave_index < static_cast<Index>(slave_surfaces->size());
         ++slave_index) {
#ifdef _OPENMP
        if (assembly_failed.load(std::memory_order_relaxed)) {
            continue;
        }
        const int thread_id = omp_get_thread_num();
#else
        const int thread_id = 0;
#endif

        auto& constraints = thread_constraints[static_cast<std::size_t>(thread_id)];

        try {
            const ID slave_surface_id =
                slave_surfaces->at(static_cast<std::size_t>(slave_index));

            if (!valid_surface_id(slave_surface_id, surfaces)) {
                continue;
            }

            const auto& slave = surfaces[static_cast<std::size_t>(slave_surface_id)];

            DynamicMatrix dual;
            DynamicVector local_support;

            if (!build_dual_basis(slave, node_coords, quadrature, dual, local_support)) {
                throw std::runtime_error(
                    "CONTACT: failed to construct dual basis for slave surface " +
                    std::to_string(slave_surface_id)
                );
            }

            const Precision slave_area   = std::max(slave->area(node_coords), Precision(0));
            const Precision slave_length = std::sqrt(slave_area);

            for (Index local_node = 0; local_node < slave->n_nodes; ++local_node) {
                MortarConstraint& constraint = constraints[slave->nodes()[local_node]];
                constraint.support += local_support(local_node);
                constraint.characteristic_length =
                    std::max(constraint.characteristic_length, slave_length);
            }

            // -----------------------------------------------------------------
            // 3. Segment every slave subtriangle against every master surface
            // -----------------------------------------------------------------
            for (const LocalTriangle& slave_local : local_triangles(slave)) {
                std::array<Vec3, 3> slave_global;
                for (std::size_t i = 0; i < slave_global.size(); ++i) {
                    slave_global[i] = slave->local_to_global(slave_local[i], node_coords);
                }

                const Vec2 slave_center =
                    (slave_local[0] + slave_local[1] + slave_local[2]) / Precision(3);

                const Vec3 plane_origin =
                    slave->local_to_global(slave_center, node_coords);

                Vec3 plane_normal = slave->normal(node_coords, slave_center);

                if (!plane_origin.allFinite() || !plane_normal.allFinite() ||
                    plane_normal.norm() <= geometry_tolerance) {
                    continue;
                }
                plane_normal.normalize();

                Vec3 first  = slave_global[1] - slave_global[0];
                Vec3 second = slave_global[2] - slave_global[0];

                first  -= first.dot(plane_normal)  * plane_normal;
                second -= second.dot(plane_normal) * plane_normal;

                Vec3 plane_tangent =
                    first.squaredNorm() >= second.squaredNorm() ? first : second;

                if (!plane_tangent.allFinite() ||
                    plane_tangent.norm() <= geometry_tolerance) {
                    continue;
                }
                plane_tangent.normalize();

                Vec3 plane_binormal = plane_normal.cross(plane_tangent);
                if (!plane_binormal.allFinite() ||
                    plane_binormal.norm() <= geometry_tolerance) {
                    continue;
                }
                plane_binormal.normalize();

                auto project = [&](const Vec3& point) -> Vec2 {
                    const Vec3 delta = point - plane_origin;
                    return Vec2(delta.dot(plane_tangent), delta.dot(plane_binormal));
                };

                PlaneTriangle slave_plane {
                    project(slave_global[0]),
                    project(slave_global[1]),
                    project(slave_global[2])
                };

                if (std::abs(cross_2d(
                        slave_plane[1] - slave_plane[0],
                        slave_plane[2] - slave_plane[0]
                    )) <= area_tolerance) {
                    continue;
                }

                SurfacePolygon<3> slave_polygon {
                    slave_plane[0],
                    slave_plane[1],
                    slave_plane[2]
                };
                if (!slave_polygon.is_ccw()) {
                    slave_polygon.flip();
                }

                std::vector<ProjectedSegment> segments;

                for (ID master_surface_id : *master_surfaces) {
                    if (!valid_surface_id(master_surface_id, surfaces)) {
                        continue;
                    }

                    const auto& master = surfaces[static_cast<std::size_t>(master_surface_id)];
                    if (surfaces_share_node(*slave, *master)) {
                        continue;
                    }

                    for (const LocalTriangle& master_local : local_triangles(master)) {
                        const Vec2 master_center =
                            (master_local[0] + master_local[1] + master_local[2]) / Precision(3);

                        Vec3 master_center_normal = master->normal(node_coords, master_center);
                        if (flip_normal) {
                            master_center_normal = -master_center_normal;
                        }

                        if (!master_center_normal.allFinite() ||
                            master_center_normal.norm() <= geometry_tolerance) {
                            continue;
                        }
                        master_center_normal.normalize();

                        if (plane_normal.dot(master_center_normal) >
                            maximum_opposing_normal_dot) {
                            continue;
                        }

                        PlaneTriangle master_plane;
                        for (std::size_t i = 0; i < master_plane.size(); ++i) {
                            const Vec3 master_global =
                                master->local_to_global(master_local[i], node_coords);
                            master_plane[i] = project(master_global);
                        }

                        if (std::abs(cross_2d(
                                master_plane[1] - master_plane[0],
                                master_plane[2] - master_plane[0]
                            )) <= area_tolerance) {
                            continue;
                        }

                        SurfacePolygon<3> master_polygon {
                            master_plane[0],
                            master_plane[1],
                            master_plane[2]
                        };
                        if (!master_polygon.is_ccw()) {
                            master_polygon.flip();
                        }

                        const auto overlap = master_polygon.intersection(slave_polygon);
                        if (overlap.size() < 3 || overlap.area() <= area_tolerance) {
                            continue;
                        }

                        ProjectedSegment segment;
                        segment.master_surface_id = master_surface_id;
                        segment.master_local      = master_local;
                        segment.master_plane      = master_plane;
                        segment.overlap           = overlap;
                        segments.push_back(std::move(segment));
                    }
                }

                if (segments.empty()) {
                    continue;
                }

                // -------------------------------------------------------------
                // 4. Integrate overlap and choose one visible master layer
                // -------------------------------------------------------------
                for (std::size_t segment_id = 0; segment_id < segments.size(); ++segment_id) {
                    const ProjectedSegment& segment = segments[segment_id];
                    const Vec2 origin = segment.overlap[0];

                    for (std::size_t fan = 1; fan + 1 < segment.overlap.size(); ++fan) {
                        const Vec2 edge_r = segment.overlap[fan]     - origin;
                        const Vec2 edge_s = segment.overlap[fan + 1] - origin;

                        const Precision triangle_jacobian =
                            std::abs(cross_2d(edge_r, edge_s));

                        if (triangle_jacobian <= area_tolerance) {
                            continue;
                        }

                        for (Index q = 0; q < quadrature.count(); ++q) {
                            const auto qp = quadrature.get_point(q);

                            const Vec2 plane_point =
                                origin + qp.r * edge_r + qp.s * edge_s;
                            const Precision weight = qp.w * triangle_jacobian;

                            const Vec3 slave_lambda = barycentric(plane_point, slave_plane);
                            if (!barycentric_inside(slave_lambda)) {
                                continue;
                            }

                            const Vec2 slave_local_qp =
                                interpolate_local(slave_local, slave_lambda);
                            const Vec3 slave_global_qp =
                                slave->local_to_global(slave_local_qp, node_coords);

                            Vec3 slave_normal_qp = slave->normal(node_coords, slave_local_qp);
                            if (!slave_global_qp.allFinite() ||
                                !slave_normal_qp.allFinite() ||
                                slave_normal_qp.norm() <= geometry_tolerance) {
                                continue;
                            }
                            slave_normal_qp.normalize();

                            std::size_t selected_segment = segments.size();
                            Precision selected_distance =
                                std::numeric_limits<Precision>::max();

                            Vec2 selected_master_local  = Vec2::Zero();
                            Vec3 selected_master_normal = Vec3::Zero();
                            Precision selected_gap      = Precision(0);

                            for (std::size_t candidate_id = 0;
                                 candidate_id < segments.size();
                                 ++candidate_id) {
                                const ProjectedSegment& candidate = segments[candidate_id];
                                const Vec3 master_lambda =
                                    barycentric(plane_point, candidate.master_plane);

                                if (!barycentric_inside(master_lambda) ||
                                    !valid_surface_id(candidate.master_surface_id, surfaces)) {
                                    continue;
                                }

                                const auto& master = surfaces[static_cast<std::size_t>(
                                    candidate.master_surface_id
                                )];

                                const Vec2 master_local_qp =
                                    interpolate_local(candidate.master_local, master_lambda);
                                const Vec3 master_global_qp =
                                    master->local_to_global(master_local_qp, node_coords);

                                Vec3 master_normal_qp =
                                    master->normal(node_coords, master_local_qp);
                                if (flip_normal) {
                                    master_normal_qp = -master_normal_qp;
                                }

                                if (!master_global_qp.allFinite() ||
                                    !master_normal_qp.allFinite() ||
                                    master_normal_qp.norm() <= geometry_tolerance) {
                                    continue;
                                }
                                master_normal_qp.normalize();

                                if (slave_normal_qp.dot(master_normal_qp) >
                                    maximum_opposing_normal_dot) {
                                    continue;
                                }

                                const Vec3 difference = slave_global_qp - master_global_qp;
                                const Precision candidate_distance = difference.norm();
                                const Precision candidate_gap =
                                    difference.dot(master_normal_qp) - clearance;

                                if (!std::isfinite(candidate_distance) ||
                                    !std::isfinite(candidate_gap)) {
                                    continue;
                                }

                                const bool better =
                                    selected_segment == segments.size() ||
                                    candidate_distance <
                                        selected_distance - projection_selection_tolerance ||
                                    (std::abs(candidate_distance - selected_distance) <=
                                         projection_selection_tolerance &&
                                     candidate_id < selected_segment);

                                if (!better) {
                                    continue;
                                }

                                selected_segment       = candidate_id;
                                selected_distance      = candidate_distance;
                                selected_master_local  = master_local_qp;
                                selected_master_normal = master_normal_qp;
                                selected_gap           = candidate_gap;
                            }

                            if (selected_segment != segment_id) {
                                continue;
                            }

                            const auto& master = surfaces[static_cast<std::size_t>(
                                segments[selected_segment].master_surface_id
                            )];

                            const DynamicVector Ns = slave->shape_function(slave_local_qp);
                            const DynamicVector Nm = master->shape_function(selected_master_local);
                            const DynamicVector Phi = dual * Ns;

                            if (!Ns.allFinite() || !Nm.allFinite() || !Phi.allFinite() ||
                                !std::isfinite(selected_gap) || !std::isfinite(weight) ||
                                weight <= Precision(0)) {
                                continue;
                            }

                            // -------------------------------------------------
                            // 5. Accumulate thread-local nodal mortar constraints
                            // -------------------------------------------------
                            for (Index constraint_node = 0;
                                 constraint_node < slave->n_nodes;
                                 ++constraint_node) {
                                const Precision phi = Phi(constraint_node);
                                if (std::abs(phi) <= Precision(1e-14)) {
                                    continue;
                                }

                                MortarConstraint& constraint =
                                    constraints[slave->nodes()[constraint_node]];

                                const Precision factor = weight * phi;

                                constraint.overlap_measure += weight * std::abs(phi);
                                constraint.gap_integral    += factor * selected_gap;
                                constraint.minimum_geometric_gap =
                                    std::min(constraint.minimum_geometric_gap, selected_gap);

                                for (Index local_node = 0;
                                     local_node < slave->n_nodes;
                                     ++local_node) {
                                    add_gradient(
                                        constraint.gradient,
                                        slave->nodes()[local_node],
                                        factor * Ns(local_node) * selected_master_normal
                                    );
                                }

                                for (Index local_node = 0;
                                     local_node < master->n_nodes;
                                     ++local_node) {
                                    add_gradient(
                                        constraint.gradient,
                                        master->nodes()[local_node],
                                        -factor * Nm(local_node) * selected_master_normal
                                    );
                                }
                            }
                        }
                    }
                }
            }
        } catch (...) {
#ifdef _OPENMP
            if (!assembly_failed.exchange(true, std::memory_order_relaxed)) {
                #pragma omp critical(fem_contact_exception)
                {
                    assembly_exception = std::current_exception();
                }
            }
#else
            throw;
#endif
        }
    }

#ifdef _OPENMP
    if (assembly_exception) {
        std::rethrow_exception(assembly_exception);
    }
#endif

    // -------------------------------------------------------------------------
    // 6. Reduce thread-local constraints
    // -------------------------------------------------------------------------
    std::unordered_map<ID, MortarConstraint> constraints;

    for (auto& local_constraints : thread_constraints) {
        for (auto& [constraint_id, local] : local_constraints) {
            MortarConstraint& constraint = constraints[constraint_id];

            constraint.support         += local.support;
            constraint.overlap_measure += local.overlap_measure;
            constraint.gap_integral    += local.gap_integral;
            constraint.characteristic_length =
                std::max(constraint.characteristic_length, local.characteristic_length);
            constraint.minimum_geometric_gap =
                std::min(constraint.minimum_geometric_gap, local.minimum_geometric_gap);

            for (const auto& [node_id, gradient] : local.gradient) {
                add_gradient(constraint.gradient, node_id, gradient);
            }
        }
    }

    // -------------------------------------------------------------------------
    // 7. Normalize constraints and assemble unilateral residual and tangent
    // -------------------------------------------------------------------------
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

        const Precision gap = constraint.gap_integral / constraint.support;
        if (!std::isfinite(gap)) {
            continue;
        }

        const auto multiplier_it = current.multipliers.find(constraint_id);
        const Precision multiplier =
            multiplier_it != current.multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        const bool geometrically_closed =
            constraint.minimum_geometric_gap < -geometry_tolerance;

        if (!geometrically_closed && !(multiplier > Precision(0))) {
            continue;
        }

        current.gaps[constraint_id] = gap;
        current.characteristic_lengths[constraint_id] =
            constraint.characteristic_length;

        const Precision pressure =
            std::max(Precision(0), multiplier - penalty * gap);

        if (!(pressure > Precision(0))) {
            continue;
        }

        for (const auto& [node_id, gradient] : constraint.gradient) {
            add_translational_force(nodal_forces, node_id, -pressure * gradient);
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
}

bool Contact::update_augmented_lagrange() const {
    State& current = state();

    for (auto it = current.multipliers.begin(); it != current.multipliers.end();) {
        if (current.gaps.find(it->first) == current.gaps.end()) {
            it = current.multipliers.erase(it);
        } else {
            ++it;
        }
    }

    bool changed = false;

    for (const auto& [constraint_id, gap] : current.gaps) {
        if (!std::isfinite(gap)) {
            continue;
        }

        const auto length_it = current.characteristic_lengths.find(constraint_id);
        const Precision characteristic_length =
            length_it != current.characteristic_lengths.end()
                ? std::max(length_it->second, Precision(0))
                : Precision(0);

        const Precision gap_tolerance =
            std::max(
                augmentation_gap_absolute_tolerance,
                augmentation_gap_relative_tolerance * characteristic_length
            );

        const auto multiplier_it = current.multipliers.find(constraint_id);
        const Precision old_multiplier =
            multiplier_it != current.multipliers.end()
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

        changed = changed || multiplier_change > multiplier_tolerance;

        if (new_multiplier > Precision(0)) {
            current.multipliers[constraint_id] = new_multiplier;
        } else if (multiplier_it != current.multipliers.end()) {
            current.multipliers.erase(constraint_id);
        }
    }

    return changed;
}

} // namespace fem::constraint
