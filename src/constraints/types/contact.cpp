/**
 * @file contact.cpp
 * @brief Implements frictionless dual-mortar surface-to-surface contact.
 *
 * The complete mortar formulation is implemented in this translation unit.
 * Every current slave subtriangle defines a slave-centered tangent plane. The
 * slave triangle and all admissible master subtriangles are projected onto that
 * plane, clipped in two dimensions and integrated over the resulting physical
 * overlap polygons. All master surfaces are tested directly; no BVH, search
 * radius or persistent master-facet ownership participates in the formulation.
 *
 * A dual interpolation is constructed on each complete slave element from the
 * full-slave mass matrix. Overlap integration accumulates one normalized nodal
 * gap constraint and its frozen-geometry gradient on every global slave mortar
 * node. For multiplier `lambda_i`, normalized gap `g_i`, support `S_i` and
 * effective penalty `epsilon`, the implemented unilateral law is
 *
 *     p_i = max(0, lambda_i - epsilon g_i),
 *
 * with contact residual and frozen-geometry tangent
 *
 *     f_c = -sum_i p_i H_i,
 *     K_c =  sum_i epsilon / S_i H_i H_i^T.
 *
 * Slave surfaces are processed independently by OpenMP workers when available.
 * Each worker owns a private mortar-constraint map which is reduced only after
 * geometric integration, so the quadrature-point hot path requires no locks.
 * Augmented-Lagrange multipliers are the only persistent physical history;
 * current gaps are transactional evaluation data used by the post-Newton
 * augmentation. The effective penalty persists across accepted increments and
 * is increased only when penetration stagnates between converged Newton solves.
 *
 * @see Contact
 * @see model::SurfaceInterface
 * @see model::SurfacePolygon
 * @see loadcase::tools::NonlinearStateManager
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

// Geometric tolerances operate in the projected contact construction. They are
// intentionally independent of the nonlinear AL convergence tolerances below.
constexpr Precision geometry_tolerance = Precision(1e-12);
constexpr Precision area_tolerance     = Precision(1e-16);
constexpr Precision tangent_tolerance  = Precision(1e-14);

// Facing solid surfaces have opposite outward normals. The weak threshold keeps
// moderately curved opposing facets while rejecting side and back faces.
constexpr Precision maximum_opposing_normal_dot    = Precision(-0.1);
constexpr Precision projection_selection_tolerance = Precision(1e-10);

// Post-Newton augmented-Lagrange convergence tolerances. The gap tolerance uses
// a characteristic slave length, while the multiplier tolerance is scaled by
// the current multiplier and the penalty force associated with that gap scale.
constexpr Precision augmentation_gap_relative_tolerance        = Precision(1e-4);
constexpr Precision augmentation_gap_absolute_tolerance        = Precision(1e-10);
constexpr Precision augmentation_multiplier_relative_tolerance = Precision(1e-6);
constexpr Precision augmentation_multiplier_absolute_tolerance = Precision(1e-10);

using LocalTriangle = std::array<Vec2, 3>;
using PlaneTriangle = std::array<Vec2, 3>;

/**
 * @brief One master subtriangle represented on the current slave tangent plane.
 *
 * Natural master coordinates are retained together with their projected plane
 * coordinates so a common-plane quadrature point can be mapped back to the
 * physical master surface. `overlap` is the clipped polygon shared with the
 * active slave subtriangle.
 */
struct ProjectedSegment {
    ID                master_surface_id = ID(-1);
    LocalTriangle     master_local {};
    PlaneTriangle     master_plane {};
    SurfacePolygon<6> overlap;
};

/**
 * @brief Integrated data of one global slave-node mortar constraint.
 *
 * `support` is the complete slave dual support `S_i`. `gap_integral` stores the
 * overlap integral `c_i = int Phi_i g_n dA`, and `gradient` stores the assembled
 * frozen-geometry vector `H_i`. `minimum_geometric_gap` is a physical
 * quadrature-point penetration gate that prevents sign-changing dual functions
 * from activating a completely open interface.
 */
struct MortarConstraint {
    Precision support               = Precision(0);
    Precision overlap_measure       = Precision(0);
    Precision gap_integral          = Precision(0);
    Precision characteristic_length = Precision(0);
    Precision minimum_geometric_gap = std::numeric_limits<Precision>::max();

    std::unordered_map<ID, Vec3> gradient;
};

/**
 * Checks whether a surface id addresses an existing surface object.
 *
 * Contact regions store ids rather than owning the corresponding surfaces. The
 * check therefore guards all later geometry access against sparse or invalid
 * entries in `ModelData::surfaces`.
 *
 * @param surface_id Candidate global surface id.
 * @param surfaces Global surface repository.
 * @return `true` when the id is in range and references a non-null surface.
 */
bool valid_surface_id(
    ID                                                surface_id,
    const std::vector<model::SurfaceInterface::Ptr>& surfaces
) {
    return surface_id >= 0 &&
           static_cast<std::size_t>(surface_id) < surfaces.size() &&
           static_cast<bool>(surfaces[static_cast<std::size_t>(surface_id)]);
}

/**
 * Tests whether two surface facets share at least one global node.
 *
 * Shared-node slave/master pairs are excluded from contact segmentation to avoid
 * self-contact contributions between coincident facets of the same local mesh.
 *
 * @param first First surface facet.
 * @param second Second surface facet.
 * @return `true` when both facets reference the same global node id.
 */
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

/**
 * Evaluates the signed two-dimensional cross product.
 *
 * The scalar result is used for projected triangle orientation, area and
 * barycentric-coordinate calculations in the common plane.
 */
Precision cross_2d(const Vec2& a, const Vec2& b) {
    return a(0) * b(1) - a(1) * b(0);
}

/**
 * Computes barycentric coordinates of a point in a projected triangle.
 *
 * All quantities are expressed in the active two-dimensional common plane. A
 * degenerate triangle returns NaNs instead of dividing by a near-zero projected
 * area; callers reject such coordinates through `barycentric_inside()`.
 *
 * @param point Point in common-plane coordinates.
 * @param triangle Projected triangle in the same coordinate system.
 * @return Barycentric weights `(lambda_0, lambda_1, lambda_2)`.
 */
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

/**
 * Checks whether barycentric coordinates represent a point inside a triangle.
 *
 * A small tolerance admits points on numerically clipped polygon boundaries.
 * Non-finite coordinates are always rejected.
 */
bool barycentric_inside(const Vec3& lambda) {
    constexpr Precision tolerance = Precision(1e-10);

    return lambda.allFinite() &&
           lambda.minCoeff() >= -tolerance &&
           lambda.maxCoeff() <= Precision(1) + tolerance;
}

/**
 * Interpolates a surface natural coordinate from triangle barycentric weights.
 *
 * @param triangle Natural coordinates of the subtriangle vertices.
 * @param lambda Barycentric weights in the corresponding projected triangle.
 * @return Natural two-dimensional surface coordinate.
 */
Vec2 interpolate_local(const LocalTriangle& triangle, const Vec3& lambda) {
    return lambda(0) * triangle[0] +
           lambda(1) * triangle[1] +
           lambda(2) * triangle[2];
}

/**
 * Triangulates the natural domain of a supported surface facet.
 *
 * Linear surfaces use a fan of their natural-domain polygon. Quadratic S6 and
 * S8 surfaces use explicit subtriangulations that include midside locations so
 * curved geometry is sampled consistently by the common-plane segmentation.
 *
 * @param surface Surface whose natural domain is subdivided.
 * @return Natural-coordinate subtriangles covering the surface domain.
 */
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

/**
 * Constructs the local dual mortar interpolation on a complete slave surface.
 *
 * The full-slave mass matrix
 *
 *     M = int_Gamma_s N_s N_s^T dA
 *
 * is integrated in the current configuration and factorized by LDLT. For linear
 * S3/S4 surfaces the nodal support is `D_i = int N_i dA = (M 1)_i`, yielding the
 * dual interpolation `Phi = D M^-1 N_s`. Quadratic S6/S8 currently retain the
 * established positive diagonal support normalized to the complete surface area.
 *
 * Degenerate, non-finite or singular mass matrices are rejected before the dual
 * basis is used by overlap integration.
 *
 * @param surface Complete slave surface facet.
 * @param node_coords Current global nodal positions.
 * @param quadrature Triangle quadrature used for full-slave integration.
 * @param dual Output matrix multiplying the primal slave shape vector.
 * @param support Output complete-slave nodal support `S_i`.
 * @return `true` when a finite dual basis could be constructed.
 */
bool build_dual_basis(
    const model::SurfaceInterface::Ptr& surface,
    const model::Field&                 node_coords,
    const math::quadrature::Quadrature& quadrature,
    DynamicMatrix&                      dual,
    DynamicVector&                      support
) {
    const Index n = surface->n_nodes;
    DynamicMatrix mass = DynamicMatrix::Zero(n, n);

    // Integrate the full-slave mortar mass matrix in physical surface measure.
    surface->integrate_triangular(
        node_coords,
        surface->local_domain_polygon(),
        quadrature,
        [&](const Vec2& local, const Vec3&, Precision weight) {
            const DynamicVector shape = surface->shape_function(local);
            mass.noalias() += weight * shape * shape.transpose();
        }
    );

    // Use the largest diagonal entry as a local scale for degeneracy checks.
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

    // Linear S3/S4 surfaces use D_i = int N_i dA. Quadratic S6/S8 surfaces
    // retain the existing positive full-area normalization until a dedicated
    // higher-order dual multiplier space is introduced.
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

    // Scale each row of M^-1 by its nodal support so dual * N_s evaluates Phi.
    dual = inverse;

    for (Index i = 0; i < n; ++i) {
        if (!(support(i) > Precision(1e-14) * scale) || !std::isfinite(support(i))) {
            return false;
        }
        dual.row(i) *= support(i);
    }

    return dual.allFinite();
}

/**
 * Adds one translational nodal force contribution to a six-DOF nodal field.
 *
 * Rotational entries are intentionally untouched because frictionless mortar
 * contact acts only through translational normal forces.
 */
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

/**
 * Accumulates one node block of a mortar constraint gradient.
 *
 * Multiple slave/master interpolation contributions may address the same global
 * node, so insertion and accumulation are performed in one place.
 */
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
 * Scatters one 3x3 translational contact tangent block into sparse triplets.
 *
 * Inactive global DOFs are skipped according to `system_nodal_dofs`. Numerically
 * negligible entries are omitted to avoid filling the sparse matrix with roundoff
 * noise from outer products of constraint gradients.
 */
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

/**
 * Constructs one surface-to-surface mortar contact definition.
 *
 * The supplied penalty is stored both as the immutable analysis starting value
 * and as the mutable effective penalty used by the current AL continuation. No
 * geometric preprocessing is performed here because overlap is reconstructed
 * from the current configuration during every assembly.
 */
Contact::Contact(
    model::SurfaceRegion::Ptr master,
    model::SurfaceRegion::Ptr slave,
    Precision                 penalty_stiffness,
    Precision                 contact_clearance,
    bool                      flip_master_normal
)
    : master_surfaces(std::move(master)),
      slave_surfaces(std::move(slave)),
      initial_penalty(penalty_stiffness),
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

/**
 * Returns the innermost active transactional contact state.
 *
 * Outside any trial the committed state is returned. The method is `const`
 * because trial/history storage is numerical solver state rather than part of
 * the immutable contact definition.
 */
Contact::State& Contact::state() const {
    return trial_states.empty() ? committed_state : trial_states.back();
}

/**
 * Opens one nested transactional contact trial.
 *
 * The current multiplier and accepted evaluation state are copied. Contact
 * geometry is not copied or frozen; each subsequent assembly reconstructs it
 * directly from the nodal positions supplied through `ModelData`.
 */
void Contact::begin_trial() const {
    trial_states.push_back(state());
}

/**
 * Commits the innermost contact trial into its parent transaction.
 *
 * A top-level commit promotes the accepted multiplier/evaluation state to the
 * persistent committed state. Nested commits replace the parent trial, which is
 * used for accepted line-search candidates and post-Newton AL refreshes.
 */
void Contact::commit_trial() const {
    logging::error(!trial_states.empty(),
        "CONTACT: no trial state to commit");

    State accepted = std::move(trial_states.back());
    trial_states.pop_back();

    if (trial_states.empty()) {
        committed_state = std::move(accepted);
    } else {
        trial_states.back() = std::move(accepted);
    }
}

/**
 * Discards the innermost transactional contact trial.
 *
 * The parent or committed state remains unchanged. Because geometry has no
 * persistent state, rollback consists only of dropping the trial multiplier and
 * accepted gap data.
 */
void Contact::rollback_trial() const {
    logging::error(!trial_states.empty(),
        "CONTACT: no trial state to roll back");

    trial_states.pop_back();
}

/**
 * Resets the numerical AL penalty state for a new nonlinear analysis.
 *
 * The effective penalty returns to the user-provided `initial_penalty`, and the
 * penetration history used to detect stagnation is cleared. This operation is
 * intentionally performed once per nonlinear analysis rather than once per load
 * increment, so a useful adapted penalty persists along the accepted load path.
 */
void Contact::reset_penalty_adaptation() const {
    penalty                           = initial_penalty;
    previous_augmentation_penetration = Precision(-1);
    previous_augmentation_changed     = false;
}

/**
 * Adapts the effective contact penalty after a converged inner Newton solve.
 *
 * The maximum represented penetration
 *
 *     d_k = max_i max(0, -g_i)
 *
 * is compared with the penetration stored after the previous successful
 * multiplier augmentation. If that augmentation changed at least one multiplier
 * but reduced penetration by less than 20 %, i.e.
 *
 *     d_k >= 0.8 d_(k-1),
 *
 * the effective penalty is increased by one decade. The value is capped at
 * `1e6 * initial_penalty`. The comparison is evaluated before the next multiplier
 * update so that update immediately uses the adapted penalty.
 *
 * @return `true` only when the effective penalty was increased.
 */
bool Contact::adapt_penalty() const {
    State& current = state();

    Precision maximum_penetration = Precision(0);
    for (const auto& [constraint_id, gap] : current.gaps) {
        (void) constraint_id;
        if (std::isfinite(gap)) {
            maximum_penetration =
                std::max(maximum_penetration, std::max(Precision(0), -gap));
        }
    }

    // Establish a new reference penetration when no preceding augmentation
    // changed the multiplier state or the interface is currently penetration-free.
    if (!previous_augmentation_changed || maximum_penetration <= Precision(0)) {
        previous_augmentation_penetration = maximum_penetration;
        return false;
    }

    const Precision previous = previous_augmentation_penetration;
    previous_augmentation_penetration = maximum_penetration;

    if (previous <= Precision(0) ||
        maximum_penetration < Precision(0.8) * previous) {
        return false;
    }

    const Precision maximum_penalty = initial_penalty * Precision(1e6);
    const Precision adapted_penalty =
        std::min(maximum_penalty, penalty * Precision(10));

    if (!(adapted_penalty > penalty) || !std::isfinite(adapted_penalty)) {
        return false;
    }

    penalty = adapted_penalty;
    return true;
}

/**
 * Records whether the most recent AL multiplier update changed contact history.
 *
 * The flag controls whether the next converged penetration is eligible for the
 * stagnation test in `adapt_penalty()`.
 */
void Contact::finish_augmentation(bool changed) const {
    previous_augmentation_changed = changed;
}

/**
 * Returns the effective penalty currently used by contact residual, tangent and
 * multiplier updates. The value may exceed the user input after stagnation-based
 * continuation.
 */
Precision Contact::current_penalty() const noexcept {
    return penalty;
}

/**
 * Assembles the current dual-mortar contact residual and frozen-geometry tangent.
 *
 * The algorithm proceeds directly from current nodal positions:
 *
 * 1. construct a complete-slave dual basis and nodal support,
 * 2. split each slave facet into natural-coordinate subtriangles,
 * 3. build one tangent plane from the current slave subtriangle,
 * 4. project and clip every admissible master subtriangle on that plane,
 * 5. integrate overlap polygons and select the nearest physical master layer at
 *    every quadrature point,
 * 6. accumulate thread-local nodal gap integrals and gradients,
 * 7. reduce worker maps and evaluate the unilateral AL law,
 * 8. scatter `f_c` and the symmetric frozen-geometry tangent `K_c`.
 *
 * The master/slave geometry and normals are recomputed every call. Derivatives of
 * segmentation, closest-layer selection and normals are deliberately omitted from
 * the tangent; only the fixed-current-geometry gap derivative is represented.
 *
 * @param system_nodal_dofs Global active-DOF ids used for tangent scattering.
 * @param model_data Model repository providing current positions and surfaces.
 * @param nodal_forces Nodal internal-force field receiving contact forces.
 * @param triplets Sparse triplet list receiving contact tangent entries.
 */
void Contact::assemble(
    SystemDofIds&     system_nodal_dofs,
    model::ModelData& model_data,
    model::NodeData&  nodal_forces,
    TripletList&      triplets
) const {
    // -------------------------------------------------------------------------
    // 1. Validate current model state and prepare thread-local accumulation
    // -------------------------------------------------------------------------
    logging::error(master_surfaces != nullptr && slave_surfaces != nullptr,
        "CONTACT: surface regions are not defined");
    logging::error(model_data.positions != nullptr,
        "CONTACT: POSITION field is not set");

    const model::Field& node_coords = *model_data.positions;
    auto&               surfaces    = model_data.surfaces;
    State&              current     = state();

    // Gaps and characteristic lengths belong to the current evaluation. The
    // multiplier map remains untouched until the post-Newton AL update.
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

    // Every worker accumulates complete mortar constraints independently. This
    // avoids locks around unordered-map and gradient updates in the nested
    // slave/master/overlap/quadrature loops.
    std::vector<std::unordered_map<ID, MortarConstraint>> thread_constraints(
        static_cast<std::size_t>(num_threads)
    );

    std::atomic_bool   assembly_failed{false};
    std::exception_ptr assembly_exception = nullptr;

    // -------------------------------------------------------------------------
    // 2. Process complete slave surfaces independently
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

            // The dual basis and support are defined on the complete slave
            // element, not separately on the later overlap segments.
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

                // The common plane is deliberately slave-centered: origin and
                // normal are evaluated at the natural center of this current
                // slave subtriangle. Master geometry does not alter the plane.
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

                // Build an orthonormal in-plane basis from the better-conditioned
                // projected slave edge. Both slave and master vertices are then
                // represented in this same two-dimensional coordinate system.
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

                // Brute-force all explicitly paired master facets. This is the
                // parameter-free correctness baseline: there is no DISTANCE gate
                // and no persistent spatial-search or partner state.
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
                // 4. Integrate overlap and select one physical master layer
                // -------------------------------------------------------------
                for (std::size_t segment_id = 0; segment_id < segments.size(); ++segment_id) {
                    const ProjectedSegment& segment = segments[segment_id];
                    const Vec2 origin = segment.overlap[0];

                    // A clipped convex polygon is integrated by a triangle fan
                    // around its first vertex. The stored quadrature is defined
                    // on the reference triangle and scaled by the 2D Jacobian.
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

                            // Recover slave natural coordinates from the current
                            // plane point and evaluate its physical position and
                            // outward normal in the deformed configuration.
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

                            // Several projected master layers can cover the same
                            // plane point. Select the physically nearest opposing
                            // layer; deterministic segment-id tie breaking avoids
                            // non-repeatable ownership at equal distance.
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

                            // Each plane point is integrated exactly once by the
                            // segment selected as the nearest physical master.
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

                                // c_i = int Phi_i g_n dA. `overlap_measure`
                                // uses |Phi_i| only as a coverage gate; it does
                                // not enter the normalized mortar equation.
                                constraint.overlap_measure += weight * std::abs(phi);
                                constraint.gap_integral    += factor * selected_gap;
                                constraint.minimum_geometric_gap =
                                    std::min(constraint.minimum_geometric_gap, selected_gap);

                                // H_i = int Phi_i B^T n dA. Slave translation
                                // enters with +Ns n and master translation with
                                // -Nm n for g_n = (x_s - x_m)^T n_m - clearance.
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
            // Exceptions may not escape an OpenMP worker. Record the first one
            // and rethrow it after the parallel region has completed.
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
    // 6. Reduce thread-local constraints into global slave-node constraints
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

        // Normalize the integrated dual gap by the complete slave support:
        // g_i = c_i / S_i.
        const Precision gap = constraint.gap_integral / constraint.support;
        if (!std::isfinite(gap)) {
            continue;
        }

        const auto multiplier_it = current.multipliers.find(constraint_id);
        const Precision multiplier =
            multiplier_it != current.multipliers.end()
                ? std::max(multiplier_it->second, Precision(0))
                : Precision(0);

        // A sign-changing dual basis can produce a negative normalized gap on a
        // physically open partial overlap. New contact therefore requires at
        // least one physically penetrating quadrature point. Existing positive
        // multiplier history may persist while the interface unloads.
        const bool geometrically_closed =
            constraint.minimum_geometric_gap < -geometry_tolerance;

        if (!geometrically_closed && !(multiplier > Precision(0))) {
            continue;
        }

        current.gaps[constraint_id] = gap;
        current.characteristic_lengths[constraint_id] =
            constraint.characteristic_length;

        // Semismooth unilateral AL law p_i = max(0, lambda_i - epsilon g_i).
        const Precision pressure =
            std::max(Precision(0), multiplier - penalty * gap);

        if (!(pressure > Precision(0))) {
            continue;
        }

        // Contact residual f_c = -sum_i p_i H_i.
        for (const auto& [node_id, gradient] : constraint.gradient) {
            add_translational_force(nodal_forces, node_id, -pressure * gradient);
        }

        // Frozen-geometry derivative:
        //   dg_i = H_i^T du / S_i,
        //   dp_i = -epsilon dg_i,
        //   K_i  = epsilon / S_i H_i H_i^T.
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

/**
 * Updates global slave-node augmented-Lagrange multipliers after Newton convergence.
 *
 * Constraints absent from the accepted current overlap lose obsolete multiplier
 * history. For every represented nodal constraint, the update
 *
 *     lambda_i <- max(0, lambda_i - epsilon g_i)
 *
 * is applied only when the gap lies outside the characteristic-length-scaled
 * tolerance. Positive multipliers may therefore unload to zero for an opening
 * interface. The return value reports whether any multiplier changed beyond the
 * corresponding numerical multiplier tolerance and is used to request another
 * Newton solve at the same load factor.
 *
 * @return `true` when at least one multiplier changed meaningfully.
 */
bool Contact::update_augmented_lagrange() const {
    State& current = state();

    // Remove history for mortar constraints that are no longer represented by
    // the current accepted overlap geometry.
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
