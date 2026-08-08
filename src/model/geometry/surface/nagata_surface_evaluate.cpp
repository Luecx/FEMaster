/**
 * @file nagata_surface_evaluate.cpp
 * @brief Implements analytical differential evaluation of G1 Nagata patches.
 *
 * Evaluation combines the homogeneous quadratic C0 patch with the rational
 * quintic G1 correction. Position, first derivatives and second derivatives are
 * computed in reduced barycentric coordinates and mapped back to the natural
 * coordinates of the originating finite-element surface.
 *
 * @see NagataSurface
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#include "nagata_surface.h"

#include "../../../core/logging.h"

#include <Eigen/Geometry>
#include <Eigen/QR>

#include <cmath>
#include <limits>

namespace fem::model {

/**
 * @brief Rebuilds and evaluates one patch position from another coordinate field.
 *
 * Contact coordinate sensitivities perturb one nodal coordinate at a time. A
 * complete `NagataSurface` reconstruction would unnecessarily rebuild all
 * patches, adjacency and the global projection BVH. This operation instead
 * reconstructs only the three positions and averaged normals controlling the
 * addressed patch, derives its local C0 and G1 coefficients, and evaluates the
 * position at the unchanged closest-point coordinates.
 *
 * Patch topology and incident source-surface entries are invariant under the
 * coordinate perturbation and are reused from this reconstruction. The method
 * therefore has work proportional to one local patch fan rather than the
 * complete selected master surface.
 *
 * @param location Valid patch and reduced barycentric coordinates.
 * @param node_coords Alternative finite nodal coordinate field.
 * @return Reconstructed position for the alternative geometry.
 */
Vec3 NagataSurface::evaluate_position(const Location& location,
                                      const Field&    node_coords) const {
    // Validate the snapshot field and addressed patch before local rebuilding.
    logging::error(valid(location),
        "NagataSurface::evaluate_position received an invalid location");
    logging::error(node_coords.domain == FieldDomain::NODE && node_coords.components >= 3,
        "NagataSurface::evaluate_position requires a nodal coordinate field");

    constexpr Precision normal_tolerance           = Precision(1e-12);
    constexpr Precision coefficient_rank_tolerance = Precision(1e-4);

    const Patch& patch = patches_[static_cast<std::size_t>(location.patch)];

    std::array<Vec3, 3> position;
    std::array<Vec3, 3> normal;
    std::array<Vec3, 3> q;
    std::array<Vec3, 3> gamma;

    // Re-evaluate the three patch vertices and their smooth averaged normals.
    // Each stored source entry identifies the exact FE nodal normal that
    // contributed during construction of the baseline surface.
    for (Index i = 0; i < 3; ++i) {
        const Vertex& vertex = vertices_[static_cast<std::size_t>(
            patch.vertices[static_cast<std::size_t>(i)])];

        logging::error(vertex.source_surfaces.size() == vertex.source_local_nodes.size()
                    && !vertex.source_surfaces.empty(),
            "NagataSurface contains inconsistent vertex-normal sources");
        logging::error(vertex.node_id >= 0
                    && static_cast<Index>(vertex.node_id) < node_coords.rows,
            "NagataSurface vertex lies outside the alternative coordinate field");

        position[static_cast<std::size_t>(i)] =
            node_coords.row_vec3(static_cast<Index>(vertex.node_id));
        normal[static_cast<std::size_t>(i)].setZero();

        logging::error(position[static_cast<std::size_t>(i)].allFinite(),
            "NagataSurface encountered an invalid perturbed vertex position");

        for (std::size_t source = 0; source < vertex.source_surfaces.size(); ++source) {
            const Index surface_index = vertex.source_surfaces[source];
            const Index local_node    = vertex.source_local_nodes[source];

            logging::error(surface_index < static_cast<Index>(source_surfaces_.size()),
                "NagataSurface contains an invalid vertex-normal source surface");

            const SurfaceInterface::Ptr& surface =
                source_surfaces_[static_cast<std::size_t>(surface_index)];
            const DynamicMatrix natural_coords = surface->node_coords_natural();
            const Vec2 local =
                natural_coords.row(static_cast<Eigen::Index>(local_node)).transpose();

            Vec3 source_normal = surface->normal(node_coords, local);

            logging::error(source_normal.allFinite()
                        && source_normal.norm() > normal_tolerance,
                "NagataSurface encountered an invalid perturbed vertex normal");

            normal[static_cast<std::size_t>(i)] += source_normal.normalized();
        }

        logging::error(normal[static_cast<std::size_t>(i)].norm() > normal_tolerance,
            "NagataSurface cannot average a perturbed vertex normal");

        normal[static_cast<std::size_t>(i)].normalize();
    }

    // Reconstruct the three quadratic edge coefficients from the perturbed
    // positions and endpoint normals.
    for (Index i = 0; i < 3; ++i) {
        const Index j = (i + 1) % 3;
        const Index k = (i + 2) % 3;
        const Vec3 edge = position[static_cast<std::size_t>(k)]
                        - position[static_cast<std::size_t>(j)];

        logging::error(edge.norm() > normal_tolerance,
            "NagataSurface perturbed patch contains a zero-length edge");

        StaticMatrix<2, 3> system;
        system.row(0) = normal[static_cast<std::size_t>(j)].transpose();
        system.row(1) = normal[static_cast<std::size_t>(k)].transpose();

        Vec2 rhs;
        rhs(0) =  normal[static_cast<std::size_t>(j)].dot(edge);
        rhs(1) = -normal[static_cast<std::size_t>(k)].dot(edge);

        Eigen::CompleteOrthogonalDecomposition<StaticMatrix<2, 3>> decomposition(system);
        decomposition.setThreshold(coefficient_rank_tolerance);

        const Vec3 curvature = decomposition.solve(rhs);

        logging::error(curvature.allFinite(),
            "NagataSurface failed to reconstruct a perturbed C0 coefficient");

        q[static_cast<std::size_t>(i)] =
            position[static_cast<std::size_t>(j)]
            + position[static_cast<std::size_t>(k)]
            - curvature;
    }

    // Reconstruct the rational G1 correction using the same Nagata edge system
    // and numerical rank convention as the complete surface constructor.
    for (Index i = 0; i < 3; ++i) {
        const Index j = (i + 1) % 3;
        const Index k = (i + 2) % 3;

        const Vec3& pj = position[static_cast<std::size_t>(j)];
        const Vec3& pk = position[static_cast<std::size_t>(k)];
        const Vec3& nj = normal[static_cast<std::size_t>(j)];
        const Vec3& nk = normal[static_cast<std::size_t>(k)];
        const Vec3& qi = q[static_cast<std::size_t>(i)];
        const Vec3& qj = q[static_cast<std::size_t>(j)];
        const Vec3& qk = q[static_cast<std::size_t>(k)];

        const Vec3 t_ij = qi - Precision(2) * pj;
        const Vec3 t_ik = qi - Precision(2) * pk;
        const Vec3 s_ij = qi - qk;
        const Vec3 s_ik = Precision(2) * pk - qj;

        const Precision c = nk.dot(t_ij);
        const Precision d = nj.dot(t_ik);

        Precision c_hat = Precision(1);
        Precision d_hat = Precision(1);

        if (std::abs(c) + std::abs(d) > normal_tolerance) {
            const Precision scale = Precision(1) / (std::abs(c) + std::abs(d));
            c_hat = scale * c;
            d_hat = scale * d;
        }

        const Precision edge_mismatch =
            c_hat * nj.dot(s_ik) + d_hat * nk.dot(s_ij);

        StaticMatrix<2, 3> system;
        system.row(0) = c_hat * nj.transpose();
        system.row(1) = d_hat * nk.transpose();

        Eigen::CompleteOrthogonalDecomposition<StaticMatrix<2, 3>> decomposition(system);
        decomposition.setThreshold(coefficient_rank_tolerance);

        gamma[static_cast<std::size_t>(i)] =
            decomposition.solve(Vec2::Constant(edge_mismatch));

        logging::error(gamma[static_cast<std::size_t>(i)].allFinite(),
            "NagataSurface failed to reconstruct a perturbed G1 coefficient");
    }

    // Evaluate the locally rebuilt position in reduced barycentric coordinates.
    const Precision b0 = Precision(1) - location.local(0) - location.local(1);
    const Precision b1 = location.local(0);
    const Precision b2 = location.local(1);

    Vec3 result =
          position[0] * b0*b0
        + position[1] * b1*b1
        + position[2] * b2*b2
        + q[0] * b1*b2
        + q[1] * b2*b0
        + q[2] * b0*b1;

    const Precision denominator = b0*b1 + b1*b2 + b2*b0;

    if (denominator > Precision(0)) {
        const Vec3 numerator =
              gamma[0] * b0*b1*b1*b2*b2
            + gamma[1] * b0*b0*b1*b2*b2
            + gamma[2] * b0*b0*b1*b1*b2;

        result += numerator / denominator;
    }

    logging::error(result.allFinite(),
        "NagataSurface produced an invalid locally reconstructed position");

    return result;
}

/**
 * @brief Evaluates one reconstructed patch location and its first two derivatives.
 *
 * With `beta_0 = 1 - xi - eta`, `beta_1 = xi` and `beta_2 = eta`, the C0 patch
 * combines vertex terms `P_i beta_i^2` and mixed terms `Q_i beta_j beta_k`.
 * The G1 contribution divides a quintic vector numerator by
 * `omega = beta_0 beta_1 + beta_1 beta_2 + beta_2 beta_0`. All quotient
 * derivatives are evaluated analytically. At a vertex, where numerator and
 * denominator both vanish, the continuous positional and first-order C0
 * extension is used.
 *
 * The normal is the normalized cross product of both reconstructed tangents.
 * Source-element coordinates are the barycentric interpolation of the stored
 * natural coordinates at the patch vertices.
 *
 * @param location Valid patch ID and reduced barycentric coordinates.
 * @return Reconstructed differential geometry and source-element mapping.
 */
NagataSurface::Evaluation NagataSurface::evaluate(const Location& location) const {
    // Validate the topological and parametric location before indexing storage.
    logging::error(location.patch < static_cast<nagata::PatchID>(patches_.size()),
        "NagataSurface::evaluate received an invalid patch ID");
    logging::error(location.local.allFinite(),
        "NagataSurface::evaluate received invalid local coordinates");

    const Patch& patch = patches_[static_cast<std::size_t>(location.patch)];

    const Precision xi  = location.local(0);
    const Precision eta = location.local(1);
    const Precision b0  = Precision(1) - xi - eta;
    const Precision b1  = xi;
    const Precision b2  = eta;

    const Precision coordinate_tolerance =
        Precision(100) * std::numeric_limits<Precision>::epsilon();

    logging::error(b0 >= -coordinate_tolerance
                && b1 >= -coordinate_tolerance
                && b2 >= -coordinate_tolerance,
        "NagataSurface::evaluate location lies outside its triangular patch");

    const Vec3& p0 = vertices_[static_cast<std::size_t>(patch.vertices[0])].position;
    const Vec3& p1 = vertices_[static_cast<std::size_t>(patch.vertices[1])].position;
    const Vec3& p2 = vertices_[static_cast<std::size_t>(patch.vertices[2])].position;

    const Vec3& q0 = patch.q[0];
    const Vec3& q1 = patch.q[1];
    const Vec3& q2 = patch.q[2];

    // ---------------------------------------------------------------------
    // Evaluate the quadratic C0 patch and its constant second derivatives
    // ---------------------------------------------------------------------

    const Vec3 position_c0 =
          p0 * b0*b0
        + p1 * b1*b1
        + p2 * b2*b2
        + q0 * b1*b2
        + q1 * b2*b0
        + q2 * b0*b1;

    const Vec3 d_xi_c0 =
        - Precision(2) * p0 * b0
        + Precision(2) * p1 * b1
        + q0 * b2
        - q1 * b2
        + q2 * (b0 - b1);

    const Vec3 d_eta_c0 =
        - Precision(2) * p0 * b0
        + Precision(2) * p2 * b2
        + q0 * b1
        + q1 * (b0 - b2)
        - q2 * b1;

    const Vec3 d2_xixi_c0 =
          Precision(2) * p0
        + Precision(2) * p1
        - Precision(2) * q2;

    const Vec3 d2_xieta_c0 =
          Precision(2) * p0
        + q0
        - q1
        - q2;

    const Vec3 d2_etaeta_c0 =
          Precision(2) * p0
        + Precision(2) * p2
        - Precision(2) * q1;

    // ---------------------------------------------------------------------
    // Evaluate the rational numerator basis through second order
    // ---------------------------------------------------------------------

    const Precision phi0 = b0*b1*b1*b2*b2;
    const Precision phi1 = b0*b0*b1*b2*b2;
    const Precision phi2 = b0*b0*b1*b1*b2;

    const Precision phi0_xi  = b1*b2*b2 * (Precision(2)*b0 - b1);
    const Precision phi0_eta = b1*b1*b2 * (Precision(2)*b0 - b2);
    const Precision phi1_xi  = b0*b2*b2 * (b0 - Precision(2)*b1);
    const Precision phi1_eta = Precision(2)*b0*b1*b2 * (b0 - b2);
    const Precision phi2_xi  = Precision(2)*b0*b1*b2 * (b0 - b1);
    const Precision phi2_eta = b0*b1*b1 * (b0 - Precision(2)*b2);

    const Precision phi0_xixi =
        Precision(2)*b2*b2 * (b0 - Precision(2)*b1);
    const Precision phi0_xieta =
        Precision(2)*b1*b2 * (Precision(2)*b0 - b1 - b2);
    const Precision phi0_etaeta =
        Precision(2)*b1*b1 * (b0 - Precision(2)*b2);

    const Precision phi1_xixi =
        Precision(2)*b2*b2 * (b1 - Precision(2)*b0);
    const Precision phi1_xieta = Precision(2)*b2
        * (b0*(b0 - b2) - b1*(Precision(2)*b0 - b2));
    const Precision phi1_etaeta = Precision(2)*b1
        * (b0*b0 + b2*b2 - Precision(4)*b0*b2);

    const Precision phi2_xixi = Precision(2)*b2
        * (b0*b0 + b1*b1 - Precision(4)*b0*b1);
    const Precision phi2_xieta = Precision(2)*b1
        * (b0*(b0 - b1) - b2*(Precision(2)*b0 - b1));
    const Precision phi2_etaeta =
        Precision(2)*b1*b1 * (b2 - Precision(2)*b0);

    const Vec3& gamma0 = patch.gamma[0];
    const Vec3& gamma1 = patch.gamma[1];
    const Vec3& gamma2 = patch.gamma[2];

    const Vec3 numerator = gamma0 * phi0 + gamma1 * phi1 + gamma2 * phi2;
    const Vec3 numerator_xi =
        gamma0 * phi0_xi + gamma1 * phi1_xi + gamma2 * phi2_xi;
    const Vec3 numerator_eta =
        gamma0 * phi0_eta + gamma1 * phi1_eta + gamma2 * phi2_eta;
    const Vec3 numerator_xixi =
        gamma0 * phi0_xixi + gamma1 * phi1_xixi + gamma2 * phi2_xixi;
    const Vec3 numerator_xieta =
        gamma0 * phi0_xieta + gamma1 * phi1_xieta + gamma2 * phi2_xieta;
    const Vec3 numerator_etaeta =
        gamma0 * phi0_etaeta + gamma1 * phi1_etaeta + gamma2 * phi2_etaeta;

    // ---------------------------------------------------------------------
    // Apply the quotient rule to the rational G1 correction
    // ---------------------------------------------------------------------

    const Precision denominator = b0*b1 + b1*b2 + b2*b0;

    Vec3 correction        = Vec3::Zero();
    Vec3 correction_xi     = Vec3::Zero();
    Vec3 correction_eta    = Vec3::Zero();
    Vec3 correction_xixi   = Vec3::Zero();
    Vec3 correction_xieta  = Vec3::Zero();
    Vec3 correction_etaeta = Vec3::Zero();

    if (denominator > Precision(0)) {
        const Precision denominator_xi  = b0 - b1;
        const Precision denominator_eta = b0 - b2;

        constexpr Precision denominator_xixi   = Precision(-2);
        constexpr Precision denominator_xieta  = Precision(-1);
        constexpr Precision denominator_etaeta = Precision(-2);

        const Precision denominator2 = denominator * denominator;
        const Precision denominator3 = denominator2 * denominator;

        correction = numerator / denominator;
        correction_xi = numerator_xi / denominator
                      - numerator * denominator_xi / denominator2;
        correction_eta = numerator_eta / denominator
                       - numerator * denominator_eta / denominator2;

        correction_xixi = numerator_xixi / denominator
            - Precision(2) * numerator_xi * denominator_xi / denominator2
            - numerator * denominator_xixi / denominator2
            + Precision(2) * numerator * denominator_xi*denominator_xi / denominator3;

        correction_xieta = numerator_xieta / denominator
            - numerator_xi * denominator_eta / denominator2
            - numerator_eta * denominator_xi / denominator2
            - numerator * denominator_xieta / denominator2
            + Precision(2) * numerator * denominator_xi*denominator_eta / denominator3;

        correction_etaeta = numerator_etaeta / denominator
            - Precision(2) * numerator_eta * denominator_eta / denominator2
            - numerator * denominator_etaeta / denominator2
            + Precision(2) * numerator * denominator_eta*denominator_eta / denominator3;
    }

    // ---------------------------------------------------------------------
    // Assemble geometry and map back to the source FE natural domain
    // ---------------------------------------------------------------------

    Evaluation result;

    result.position        = position_c0 + correction;
    result.jacobian.col(0) = d_xi_c0 + correction_xi;
    result.jacobian.col(1) = d_eta_c0 + correction_eta;
    result.d2_xixi         = d2_xixi_c0 + correction_xixi;
    result.d2_xieta        = d2_xieta_c0 + correction_xieta;
    result.d2_etaeta       = d2_etaeta_c0 + correction_etaeta;

    const Vec3 normal = result.jacobian.col(0).cross(result.jacobian.col(1));

    logging::error(normal.allFinite() && normal.norm() > Precision(1e-12),
        "NagataSurface::evaluate encountered a degenerate patch Jacobian");

    result.normal = normal.normalized();

    logging::error(patch.surface < static_cast<Index>(source_surfaces_.size()),
        "NagataSurface contains an invalid source-surface reference");

    result.surface = source_surfaces_[static_cast<std::size_t>(patch.surface)].get();
    result.element_local = patch.element_local[0] * b0
                         + patch.element_local[1] * b1
                         + patch.element_local[2] * b2;

    return result;
}

/**
 * @brief Tests whether a tracked location addresses this reconstruction.
 *
 * Patch numbering depends only on the ordered source topology and therefore
 * remains stable while the nodal geometry changes between contact assemblies.
 * The local coordinates must remain finite but may lie on a closed patch edge.
 *
 * @param location Location to validate.
 * @return True for a finite location with an existing patch identifier.
 */
bool NagataSurface::valid(const Location& location) const {
    return location.patch < static_cast<nagata::PatchID>(patches_.size())
        && location.local.allFinite();
}

/**
 * @brief Returns the connected master-surface component containing a location.
 *
 * Internal triangular charts and originating FE surfaces share one component
 * whenever patch adjacency connects them. Contact can therefore track freely
 * across those internal boundaries without treating a patch change as a
 * discrete partner switch.
 *
 * @param location Valid reconstructed surface location.
 * @return Connected-component identifier assigned during reconstruction.
 */
nagata::ComponentID NagataSurface::component(const Location& location) const {
    logging::error(valid(location),
        "NagataSurface::component received an invalid location");

    return patches_[static_cast<std::size_t>(location.patch)].component;
}

} // namespace fem::model
