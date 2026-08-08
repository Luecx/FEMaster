/**
 * @file neo_hooke_elasticity.cpp
 * @brief Implements compressible isotropic Neo-Hookean elasticity.
 *
 * The implementation evaluates the three-dimensional PK2 response from the
 * right Cauchy-Green tensor and differentiates it analytically with respect to
 * Green-Lagrange strain. Axial and shell plane-stress reductions determine the
 * missing transverse deformation by Newton iteration and consistently condense
 * the full material tangent.
 *
 * All constitutive calls are state-neutral because the Neo-Hookean law contains
 * no history variables; the state pointer is accepted only through the common
 * `Elasticity` interface.
 *
 * @see NeoHookeElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "neo_hooke_elasticity.h"

#include "../core/logging.h"
#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

#include <Eigen/LU>

#include <array>
#include <cmath>

namespace fem::material {

/**
 * Constructs the compressible Neo-Hookean material and derives its infinitesimal
 * elastic moduli.
 *
 * The potential parameters must be positive. They define `mu = 2 C10`,
 * `K = 2/D1` and `lambda = K - 2 mu/3`, which are reused by the linearized
 * response and plane-stress initial guesses.
 *
 * @param c10_in Isochoric energy coefficient.
 * @param d1_in Volumetric penalty parameter.
 */
NeoHookeElasticity::NeoHookeElasticity(Precision c10_in, Precision d1_in)
    : c10(c10_in),
      d1 (d1_in) {
    logging::error(c10 > Precision(0),
                   "NEO_HOOKE: C10 must be positive");
    logging::error(d1 > Precision(0),
                   "NEO_HOOKE: D1 must be positive");

    mu          = Precision(2) * c10;
    bulk        = Precision(2) / d1;
    lame_lambda = bulk - Precision(2) * mu / Precision(3);
}

bool NeoHookeElasticity::supports_axial_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_axial_green_lagrange() const {
    return true;
}

bool NeoHookeElasticity::supports_volume_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_volume_green_lagrange() const {
    return true;
}

bool NeoHookeElasticity::supports_shell_integration_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

/**
 * Builds the infinitesimal three-dimensional tangent of the finite-strain
 * potential at the undeformed state.
 *
 * @return Isotropic engineering-Voigt tangent defined by `lambda` and `mu`.
 */
Mat6 NeoHookeElasticity::linear_tangent() const {
    Mat6 tangent;
    const Precision C11 = lame_lambda + Precision(2) * mu;
    tangent <<
        C11         , lame_lambda , lame_lambda , Precision(0), Precision(0), Precision(0),
        lame_lambda , C11         , lame_lambda , Precision(0), Precision(0), Precision(0),
        lame_lambda , lame_lambda , C11         , Precision(0), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), mu,           Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), mu,           Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), mu;
    return tangent;
}

/**
 * Builds the infinitesimal shell plane-stress tangent implied by the
 * Neo-Hookean bulk and shear moduli.
 *
 * @return Five-component tangent ordered `[11,22,12,13,23]`.
 */
Mat5 NeoHookeElasticity::linear_shell_tangent() const {
    const Precision youngs  = Precision(9) * bulk * mu / (Precision(3) * bulk + mu);
    const Precision poisson = (Precision(3) * bulk - Precision(2) * mu)
                            / (Precision(2) * (Precision(3) * bulk + mu));
    const Precision scalar  = youngs / (Precision(1) - poisson * poisson);

    Mat5 tangent = Mat5::Zero();
    tangent(0, 0) = scalar;
    tangent(0, 1) = scalar * poisson;
    tangent(1, 0) = scalar * poisson;
    tangent(1, 1) = scalar;
    tangent(2, 2) = mu;
    tangent(3, 3) = mu;
    tangent(4, 4) = mu;
    return tangent;
}

void NeoHookeElasticity::evaluate(const AxialStrainLinearized& strain,
                                  Precision*                   state,
                                  AxialStressCauchy&           stress,
                                  Precision&                   tangent) const {
    (void) state;

    const Precision youngs = Precision(9) * bulk * mu / (Precision(3) * bulk + mu);
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

/**
 * Evaluates the finite-strain uniaxial response with traction-free lateral
 * surfaces.
 *
 * The supplied Green-Lagrange strain determines the axial right Cauchy-Green
 * component `C11`. Isotropy permits both lateral components to share one
 * logarithmic unknown. Newton iteration drives the lateral PK2 stress to zero,
 * after which the axial tangent is condensed from the full three-dimensional
 * material operator.
 *
 * The logarithmic lateral variable keeps the transverse stretch positive. A
 * non-positive axial metric or a failed local plane-stress solve is reported as
 * a material error. The Neo-Hookean model is stateless, so `state` is not
 * modified.
 *
 * @param strain Axial Green-Lagrange strain in the material direction.
 * @param state Unused material-point state row retained by the common interface.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Consistent derivative of axial PK2 stress with respect to the
 *                axial Green-Lagrange strain.
 */
void NeoHookeElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                  Precision*                      state,
                                  AxialStressPK2&                 stress,
                                  Precision&                      tangent) const {
    (void) state;

    const Precision c       = Precision(1) + Precision(2) * strain.value();
    const Precision poisson = (Precision(3) * bulk - Precision(2) * mu)
                            / (Precision(2) * (Precision(3) * bulk + mu));

    logging::error(c > Precision(0),
                   "NEO_HOOKE: non-positive axial right Cauchy-Green component");

    Precision log_x                  = -poisson * std::log(c);
    bool      plane_stress_converged = false;

    Mat3 full_stress;
    Mat6 full_tangent;

    // Solve the identical lateral metric components from the zero transverse-
    // stress condition while preserving their positivity in logarithmic form
    for (Index i = 0; i < 30; ++i) {
        const Precision x = std::exp(log_x);
        Mat3 C            = Mat3::Zero();
        C(0, 0)           = c;
        C(1, 1)           = x;
        C(2, 2)           = x;

        evaluate_full(C, full_stress, full_tangent);

        const Precision residual = full_stress(1, 1);
        const Precision slope    = Precision(0.5) * x
                                 * (full_tangent(1, 1) + full_tangent(1, 2));
        const Precision delta    = -residual / slope;

        log_x += delta;

        if (std::abs(delta) < Precision(1e-12)) {
            plane_stress_converged = true;
            break;
        }
    }

    // Re-evaluate the converged three-dimensional state for output and
    // consistent tangent condensation
    const Precision x = std::exp(log_x);
    Mat3 C            = Mat3::Zero();
    C(0, 0)           = c;
    C(1, 1)           = x;
    C(2, 2)           = x;
    evaluate_full(C, full_stress, full_tangent);

    logging::error(plane_stress_converged &&
                   std::abs(full_stress(1, 1)) <= Precision(1e-10) * (mu + bulk),
                   "NEO_HOOKE: axial plane-stress iteration did not converge");

    // Condense both traction-free lateral components by their Schur complement
    const Mat2 lateral_tangent = full_tangent.template block<2, 2>(1, 1);
    const Vec2 axial_to_lateral(full_tangent(1, 0), full_tangent(2, 0));
    const Vec2 lateral_to_axial(full_tangent(0, 1), full_tangent(0, 2));

    stress.value() = full_stress(0, 0);
    tangent        = full_tangent(0, 0)
                   - (lateral_to_axial.transpose()
                   * lateral_tangent.inverse()
                   * axial_to_lateral)(0, 0);
}

void NeoHookeElasticity::evaluate(const VolumeStrainLinearized& strain,
                                  Precision*                    state,
                                  VolumeStressCauchy&           stress,
                                  Mat6&                         tangent) const {
    (void) state;

    tangent        = linear_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates the full three-dimensional finite-strain material response.
 *
 * Green-Lagrange strain is converted to the right Cauchy-Green tensor through
 * `C = I + 2 E`. The full constitutive kernel returns PK2 stress and the
 * consistent material tangent in engineering-Voigt ordering. No dimensional
 * reduction or stress transformation is applied here.
 *
 * @param strain Green-Lagrange strain in the material basis.
 * @param state Unused material-point state row retained by the common interface.
 * @param stress Second Piola-Kirchhoff stress in the material basis.
 * @param tangent Consistent derivative `dS/dE`.
 */
void NeoHookeElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                  Precision*                       state,
                                  VolumeStressPK2&                 stress,
                                  Mat6&                            tangent) const {
    (void) state;

    const Mat3 C = Mat3::Identity() + Precision(2) * strain.tensor();
    Mat3       full_stress;

    evaluate_full(C, full_stress, tangent);
    stress = VolumeStressPK2(full_stress);
}

void NeoHookeElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                  Precision*                           state,
                                  ShellMaterialStressCauchy&            stress,
                                  Mat5&                                 tangent) const {
    (void) state;

    tangent         = linear_shell_tangent();
    stress.values() = tangent * strain.values();
}

/**
 * Evaluates the finite-strain shell material response under plane stress.
 *
 * The five supplied shell strain components define the in-plane and transverse-
 * shear entries of the right Cauchy-Green tensor. The missing thickness metric
 * is written as the minimum positive-definite value plus a positive logarithmic
 * unknown. Newton iteration enforces `S33 = 0` without violating positive
 * definiteness.
 *
 * After convergence, the shell PK2 components are extracted in the ordering
 * `[S11,S22,S12,S13,S23]`. The corresponding five-by-five tangent is the Schur
 * complement of the eliminated thickness component and is therefore consistent
 * with the local plane-stress solve. The state row remains unchanged.
 *
 * @param strain Green-Lagrange shell material strain.
 * @param state Unused material-point state row retained by the common interface.
 * @param stress Plane-stress shell PK2 components.
 * @param tangent Condensed consistent shell material tangent.
 */
void NeoHookeElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                  Precision*                              state,
                                  ShellMaterialStressPK2&                 stress,
                                  Mat5&                                   tangent) const {
    (void) state;

    Mat3 C = Mat3::Identity();
    C(0, 0) = Precision(1) + Precision(2) * strain.values()(0);
    C(1, 1) = Precision(1) + Precision(2) * strain.values()(1);
    C(0, 1) = strain.values()(2);
    C(1, 0) = strain.values()(2);
    C(0, 2) = strain.values()(3);
    C(2, 0) = strain.values()(3);
    C(1, 2) = strain.values()(4);
    C(2, 1) = strain.values()(4);

    const Mat2 in_plane = C.template block<2, 2>(0, 0);
    const Vec2 shear_column(C(0, 2), C(1, 2));
    const Precision det_in_plane = in_plane.determinant();

    logging::error(det_in_plane > Precision(0),
                   "NEO_HOOKE: non-positive shell in-plane right Cauchy-Green determinant");

    // The Schur complement of the in-plane block gives the smallest admissible
    // C33 for a positive-definite right Cauchy-Green tensor
    const Precision min_c33 = (shear_column.transpose() * in_plane.inverse() * shear_column)(0, 0);
    const Precision poisson = (Precision(3) * bulk - Precision(2) * mu)
                            / (Precision(2) * (Precision(3) * bulk + mu));

    Precision log_d                  = -poisson / (Precision(1) - poisson) * std::log(det_in_plane);
    bool      plane_stress_converged = false;

    Mat3 full_stress;
    Mat6 full_tangent;

    // Solve the traction-free thickness condition in a positive logarithmic
    // distance from the admissibility boundary
    for (Index i = 0; i < 30; ++i) {
        const Precision d = std::exp(log_d);
        C(2, 2)           = min_c33 + d;
        evaluate_full(C, full_stress, full_tangent);

        const Precision residual = full_stress(2, 2);
        const Precision slope    = Precision(0.5) * d * full_tangent(2, 2);
        const Precision delta    = -residual / slope;

        log_d += delta;

        if (std::abs(delta) < Precision(1e-12)) {
            plane_stress_converged = true;
            break;
        }
    }

    // Re-evaluate the converged three-dimensional response before extracting
    // shell stresses and condensing the material tangent
    const Precision d   = std::exp(log_d);
    const Precision c33 = min_c33 + d;

    C(2, 2) = c33;
    evaluate_full(C, full_stress, full_tangent);

    logging::error(plane_stress_converged &&
                   std::abs(full_stress(2, 2)) <= Precision(1e-10) * (mu + bulk),
                   "NEO_HOOKE: shell plane-stress thickness iteration did not converge");

    stress.values() << full_stress(0, 0),
                       full_stress(1, 1),
                       full_stress(0, 1),
                       full_stress(0, 2),
                       full_stress(1, 2);

    // Map the full engineering-Voigt components to shell ordering and eliminate
    // the thickness normal component by a Schur complement
    const std::array<Index, 5> shell_rows { 0, 1, 5, 4, 3 };
    const std::array<Index, 5> shell_cols { 0, 1, 5, 4, 3 };

    tangent.setZero();
    for (Index row = 0; row < 5; ++row) {
        for (Index col = 0; col < 5; ++col) {
            tangent(row, col) =
                full_tangent(shell_rows[row], shell_cols[col])
              - full_tangent(shell_rows[row], 2)
              * full_tangent(2, shell_cols[col])
              / full_tangent(2, 2);
        }
    }
}

/**
 * Evaluates the full compressible Neo-Hookean PK2 response and tangent.
 *
 * The stress follows from the isochoric-volumetric split of
 *
 *     W = C10 (J^(-2/3) I1 - 3) + (J - 1)^2 / D1.
 *
 * The tangent is assembled column-wise by applying one unit Green-Lagrange
 * variation in engineering-Voigt ordering. Shear basis variations contain
 * one-half in each symmetric tensor position so the resulting columns are the
 * derivatives with respect to engineering shear strain. The input tensor must
 * be positive definite; a non-positive determinant is rejected.
 *
 * @param C Right Cauchy-Green tensor in the material basis.
 * @param stress Second Piola-Kirchhoff stress tensor.
 * @param tangent Consistent material tangent `dS/dE` in engineering-Voigt form.
 */
void NeoHookeElasticity::evaluate_full(const Mat3& C, Mat3& stress, Mat6& tangent) const {
    const Precision det_c = C.determinant();

    logging::error(det_c > Precision(0),
                   "NEO_HOOKE: non-positive right Cauchy-Green determinant");

    const Precision J                = std::sqrt(det_c);
    const Precision first_invariant  = C.trace();
    const Precision mean_invariant   = first_invariant / Precision(3);
    const Precision deviatoric_scale = Precision(2) * c10 * std::pow(J, Precision(-2) / Precision(3));
    const Precision volumetric_scale = Precision(2) * J * (J - Precision(1)) / d1;
    const Mat3      C_inv             = C.inverse();

    // Evaluate the isochoric and volumetric PK2 stress contributions
    stress = deviatoric_scale * (Mat3::Identity() - mean_invariant * C_inv)
           + volumetric_scale * C_inv;

    // Differentiate the stress analytically for each independent
    // Green-Lagrange strain component
    tangent.setZero();
    for (Index col = 0; col < 6; ++col) {
        Mat3 dE = Mat3::Zero();

        if (col < 3) {
            dE(col, col) = Precision(1);
        } else if (col == 3) {
            dE(1, 2) = Precision(0.5);
            dE(2, 1) = Precision(0.5);
        } else if (col == 4) {
            dE(0, 2) = Precision(0.5);
            dE(2, 0) = Precision(0.5);
        } else {
            dE(0, 1) = Precision(0.5);
            dE(1, 0) = Precision(0.5);
        }

        const Precision inverse_contraction = (C_inv * dE).trace();
        const Precision strain_trace        = dE.trace();
        const Mat3      inverse_derivative  = C_inv * dE * C_inv;

        const Mat3 deviatoric_derivative = deviatoric_scale * (
            -Precision(2) / Precision(3) * inverse_contraction
                * (Mat3::Identity() - mean_invariant * C_inv)
            -Precision(2) / Precision(3) * strain_trace * C_inv
            +Precision(2) * mean_invariant * inverse_derivative
        );

        const Precision volumetric_derivative = Precision(2) * J
                                              * (Precision(2) * J - Precision(1))
                                              / d1 * inverse_contraction;
        const Mat3 dS = deviatoric_derivative
                      + volumetric_derivative * C_inv
                      - Precision(2) * volumetric_scale * inverse_derivative;

        tangent.col(col) = VolumeStressPK2(dS).voigt();
    }
}

} // namespace fem::material
