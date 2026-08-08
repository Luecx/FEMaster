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

#include <limits>

namespace fem::model {

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
