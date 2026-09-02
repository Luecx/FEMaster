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
 * Tangent output is optional through the common `Elasticity` interface. A full
 * three-dimensional stress-only query therefore skips tangent construction.
 * Axial and shell reductions still form the three-dimensional tangent internally
 * because their local Newton equations require the appropriate derivatives.
 *
 * All constitutive calls are state-neutral because the Neo-Hookean law contains
 * no history variables; the state pointers are accepted only through the common
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
 * The potential parameters must be positive. They define
 *
 *     mu     = 2 C10
 *     K      = 2 / D1
 *     lambda = K - 2 mu / 3,
 *
 * which are reused by the linearized response and plane-stress initial guesses.
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

    // Derive the infinitesimal material constants once from the energy parameters.
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
 * The normal block is formed from the Lamé constants and the engineering shear
 * entries equal `mu` because the strain vector stores engineering shear strains.
 *
 * @return Isotropic engineering-Voigt tangent defined by `lambda` and `mu`.
 */
Mat6 NeoHookeElasticity::linear_tangent() const {
    const Precision c11 = lame_lambda + Precision(2) * mu;

    Mat6 tangent;
    tangent <<
        c11,          lame_lambda, lame_lambda, Precision(0), Precision(0), Precision(0),
        lame_lambda,  c11,         lame_lambda, Precision(0), Precision(0), Precision(0),
        lame_lambda,  lame_lambda, c11,         Precision(0), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), mu,         Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), mu,         Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), mu;
    return tangent;
}

/**
 * Builds the infinitesimal shell plane-stress tangent implied by the
 * Neo-Hookean bulk and shear moduli.
 *
 * Young's modulus and Poisson's ratio are recovered from `K` and `mu`, after
 * which the in-plane normal block follows the ordinary isotropic plane-stress
 * relation. The in-plane and transverse engineering shear terms remain `mu`.
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

/**
 * Evaluates the infinitesimal axial response implied by the Neo-Hookean
 * parameters.
 *
 * The one-dimensional constitutive relation is
 *
 *     sigma = E epsilon
 *
 * with `E` recovered from the stored bulk and shear moduli. The derivative is
 * returned only when the caller provides tangent storage.
 *
 * @param strain Infinitesimal axial strain.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Axial Cauchy stress.
 * @param tangent Optional derivative `d sigma/d epsilon`.
 */
void NeoHookeElasticity::evaluate(const AxialStrainLinearized& strain,
                                  const Precision*             old_state,
                                  Precision*                   new_state,
                                  AxialStressCauchy&           stress,
                                  Precision*                   tangent) const {
    (void) old_state;
    (void) new_state;

    // Recover the infinitesimal Young's modulus from K and mu.
    const Precision youngs = Precision(9) * bulk * mu / (Precision(3) * bulk + mu);

    // Stress does not require tangent storage for this scalar linear relation.
    stress.value() = youngs * strain.value();

    if (tangent != nullptr) {
        *tangent = youngs;
    }
}

/**
 * Evaluates the finite-strain uniaxial response with traction-free lateral
 * surfaces.
 *
 * The supplied Green-Lagrange strain determines the axial right Cauchy-Green
 * component
 *
 *     C11 = 1 + 2 E11.
 *
 * Isotropy permits both lateral components to share one positive metric value
 * `x`. Its logarithm is used as the local Newton unknown so that `x > 0` is
 * preserved automatically. The local equation is
 *
 *     S22(C11, x, x) = 0.
 *
 * The full three-dimensional tangent is required by this local Newton solve even
 * when the caller does not request the final axial tangent. After convergence,
 * the optional one-dimensional derivative is obtained by eliminating both
 * traction-free lateral strains through a Schur complement.
 *
 * @param strain Axial Green-Lagrange strain in the material direction.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Optional consistent derivative `dS11/dE11`.
 */
void NeoHookeElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                  const Precision*                old_state,
                                  Precision*                      new_state,
                                  AxialStressPK2&                 stress,
                                  Precision*                      tangent) const {
    (void) old_state;
    (void) new_state;

    // Convert the prescribed axial Green-Lagrange strain into C11.
    const Precision c       = Precision(1) + Precision(2) * strain.value();
    const Precision poisson = (Precision(3) * bulk - Precision(2) * mu)
                            / (Precision(2) * (Precision(3) * bulk + mu));

    logging::error(c > Precision(0),
        "NEO_HOOKE: non-positive axial right Cauchy-Green component");

    // Use the infinitesimal Poisson response as the initial logarithmic lateral metric.
    Precision log_x                  = -poisson * std::log(c);
    bool      plane_stress_converged = false;

    Mat3 full_stress;
    Mat6 full_tangent;

    // Solve S22 = S33 = 0. Isotropy makes both equations identical, so one local
    // scalar unknown is sufficient. The derivative with respect to log(x) is
    //
    //     dS22/dlog(x)
    //         = dS22/dE22 dE22/dlog(x) + dS22/dE33 dE33/dlog(x)
    //         = 0.5 x (C2222 + C2233).
    for (Index iteration = 0; iteration < 30; ++iteration) {
        const Precision x = std::exp(log_x);

        Mat3 C  = Mat3::Zero();
        C(0, 0) = c;
        C(1, 1) = x;
        C(2, 2) = x;

        evaluate_full(C, full_stress, &full_tangent);

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

    // Re-evaluate the converged three-dimensional state. The local tangent is
    // still needed for the convergence check and for optional condensation.
    const Precision x = std::exp(log_x);

    Mat3 C  = Mat3::Zero();
    C(0, 0) = c;
    C(1, 1) = x;
    C(2, 2) = x;

    evaluate_full(C, full_stress, &full_tangent);

    logging::error(plane_stress_converged
                   && std::abs(full_stress(1, 1)) <= Precision(1e-10) * (mu + bulk),
        "NEO_HOOKE: axial plane-stress iteration did not converge");

    // Return the converged axial PK2 stress in the reference material direction.
    stress.value() = full_stress(0, 0);

    if (tangent != nullptr) {
        // Condense the two traction-free lateral components:
        //
        //     C_axial = C_aa - C_ab C_bb^-1 C_ba.
        const Mat2 lateral_tangent = full_tangent.template block<2, 2>(1, 1);
        const Vec2 axial_to_lateral(full_tangent(1, 0), full_tangent(2, 0));
        const Vec2 lateral_to_axial(full_tangent(0, 1), full_tangent(0, 2));

        *tangent = full_tangent(0, 0)
                 - (lateral_to_axial.transpose()
                 * lateral_tangent.inverse()
                 * axial_to_lateral)(0, 0);
    }
}

/**
 * Evaluates the infinitesimal three-dimensional Cauchy response.
 *
 * The constant linearized material operator is required for the stress
 * multiplication itself. It is copied to the optional tangent output only when
 * requested.
 *
 * @param strain Infinitesimal volume strain in engineering-Voigt ordering.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Cauchy stress in the material basis.
 * @param tangent Optional infinitesimal material tangent.
 */
void NeoHookeElasticity::evaluate(const VolumeStrainLinearized& strain,
                                  const Precision*              old_state,
                                  Precision*                    new_state,
                                  VolumeStressCauchy&           stress,
                                  Mat6*                         tangent) const {
    (void) old_state;
    (void) new_state;

    const Mat6 material_tangent = linear_tangent();
    stress.voigt() = material_tangent * strain.voigt();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

/**
 * Evaluates the full three-dimensional finite-strain material response.
 *
 * Green-Lagrange strain is converted to the right Cauchy-Green tensor through
 *
 *     C = I + 2 E.
 *
 * PK2 stress is always evaluated. The analytic derivative `dS/dE` is assembled
 * only when `tangent` is non-null, which makes this overload suitable for
 * residual-only nonlinear assembly.
 *
 * @param strain Green-Lagrange strain in the material basis.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Second Piola-Kirchhoff stress in the material basis.
 * @param tangent Optional consistent derivative `dS/dE`.
 */
void NeoHookeElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                  const Precision*                 old_state,
                                  Precision*                       new_state,
                                  VolumeStressPK2&                 stress,
                                  Mat6*                            tangent) const {
    (void) old_state;
    (void) new_state;

    // Convert the work-conjugate Green-Lagrange strain to its metric tensor.
    const Mat3 C = Mat3::Identity() + Precision(2) * strain.tensor();

    Mat3 full_stress;
    evaluate_full(C, full_stress, tangent);

    stress = VolumeStressPK2(full_stress);
}

/**
 * Evaluates the infinitesimal five-component shell response.
 *
 * @param strain Linearized shell material strain.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Plane-stress shell Cauchy stress.
 * @param tangent Optional reduced shell tangent.
 */
void NeoHookeElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                  const Precision*                     old_state,
                                  Precision*                           new_state,
                                  ShellMaterialStressCauchy&            stress,
                                  Mat5*                                tangent) const {
    (void) old_state;
    (void) new_state;

    const Mat5 material_tangent = linear_shell_tangent();
    stress.values() = material_tangent * strain.values();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

/**
 * Evaluates the finite-strain shell material response under plane stress.
 *
 * The five supplied shell strain components define the in-plane and transverse-
 * shear entries of the right Cauchy-Green tensor. The missing thickness metric
 * is written as the minimum positive-definite value plus a positive logarithmic
 * unknown. Newton iteration enforces
 *
 *     S33 = 0
 *
 * without violating positive definiteness.
 *
 * The local thickness solve requires the full three-dimensional tangent even for
 * a stress-only caller. After convergence, the five PK2 components are extracted
 * in shell ordering. The reduced five-by-five tangent is formed by a Schur
 * complement only when requested by the caller.
 *
 * @param strain Green-Lagrange shell material strain.
 * @param old_state Unused committed state row.
 * @param new_state Unused trial state row.
 * @param stress Plane-stress shell PK2 components.
 * @param tangent Optional condensed consistent shell tangent.
 */
void NeoHookeElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                  const Precision*                        old_state,
                                  Precision*                              new_state,
                                  ShellMaterialStressPK2&                 stress,
                                  Mat5*                                   tangent) const {
    (void) old_state;
    (void) new_state;

    // Reconstruct every prescribed entry of C = I + 2 E. Engineering shear
    // strains appear once in each symmetric off-diagonal position of C.
    Mat3 C = Mat3::Identity();
    C(0, 0) = Precision(1) + Precision(2) * strain.values()(0);
    C(1, 1) = Precision(1) + Precision(2) * strain.values()(1);
    C(0, 1) = strain.values()(2);
    C(1, 0) = strain.values()(2);
    C(0, 2) = strain.values()(3);
    C(2, 0) = strain.values()(3);
    C(1, 2) = strain.values()(4);
    C(2, 1) = strain.values()(4);

    // Positive definiteness of C requires the in-plane principal block to be
    // positive definite before a thickness metric can make the complete tensor
    // admissible.
    const Mat2      in_plane       = C.template block<2, 2>(0, 0);
    const Vec2      shear_column   (C(0, 2), C(1, 2));
    const Precision det_in_plane   = in_plane.determinant();

    logging::error(det_in_plane > Precision(0),
        "NEO_HOOKE: non-positive shell in-plane right Cauchy-Green determinant");

    // The Schur complement gives the lower admissibility bound for C33:
    //
    //     C33 > c^T Caa^-1 c
    //
    // so parameterize C33 = min_c33 + exp(log_d).
    const Precision min_c33 =
        (shear_column.transpose() * in_plane.inverse() * shear_column)(0, 0);
    const Precision poisson = (Precision(3) * bulk - Precision(2) * mu)
                            / (Precision(2) * (Precision(3) * bulk + mu));

    Precision log_d                  = -poisson / (Precision(1) - poisson)
                                     * std::log(det_in_plane);
    bool      plane_stress_converged = false;

    Mat3 full_stress;
    Mat6 full_tangent;

    // Solve the traction-free thickness condition in a positive logarithmic
    // distance from the admissibility boundary. Since
    //
    //     E33 = 0.5 (C33 - 1)
    //     dE33/dlog(d) = 0.5 d,
    //
    // the local Newton slope is 0.5 d C3333.
    for (Index iteration = 0; iteration < 30; ++iteration) {
        const Precision d = std::exp(log_d);
        C(2, 2)           = min_c33 + d;

        evaluate_full(C, full_stress, &full_tangent);

        const Precision residual = full_stress(2, 2);
        const Precision slope    = Precision(0.5) * d * full_tangent(2, 2);
        const Precision delta    = -residual / slope;

        log_d += delta;

        if (std::abs(delta) < Precision(1e-12)) {
            plane_stress_converged = true;
            break;
        }
    }

    // Re-evaluate the converged three-dimensional response. The tangent remains
    // necessary for the local convergence check and optional shell condensation.
    const Precision d   = std::exp(log_d);
    const Precision c33 = min_c33 + d;

    C(2, 2) = c33;
    evaluate_full(C, full_stress, &full_tangent);

    logging::error(plane_stress_converged
                   && std::abs(full_stress(2, 2)) <= Precision(1e-10) * (mu + bulk),
        "NEO_HOOKE: shell plane-stress thickness iteration did not converge");

    // Extract the retained shell stress components from volume ordering
    // [11,22,33,23,13,12] into [11,22,12,13,23].
    stress.values() << full_stress(0, 0),
                       full_stress(1, 1),
                       full_stress(0, 1),
                       full_stress(0, 2),
                       full_stress(1, 2);

    if (tangent != nullptr) {
        // Eliminate the thickness-normal strain by the Schur complement
        //
        //     C_shell = C_aa - C_a3 C_3a / C_33.
        const std::array<Index, 5> retained { 0, 1, 5, 4, 3 };

        tangent->setZero();
        for (Index row = 0; row < 5; ++row) {
            for (Index column = 0; column < 5; ++column) {
                (*tangent)(row, column) =
                    full_tangent(retained[row], retained[column])
                  - full_tangent(retained[row], 2)
                  * full_tangent(2, retained[column])
                  / full_tangent(2, 2);
            }
        }
    }
}

/**
 * Evaluates the full compressible Neo-Hookean PK2 response and optional tangent.
 *
 * The stress follows from the isochoric-volumetric split
 *
 *     W = C10 (J^(-2/3) I1 - 3) + (J - 1)^2 / D1.
 *
 * With
 *
 *     J = sqrt(det(C))
 *     I1 = tr(C),
 *
 * the PK2 stress is assembled directly from `C^-1`. If tangent output is
 * requested, the stress expression is differentiated analytically column by
 * column with respect to the six engineering-Voigt Green-Lagrange components.
 * Shear basis variations contain one half in each symmetric tensor position.
 *
 * @param C Right Cauchy-Green tensor in the material basis.
 * @param stress Second Piola-Kirchhoff stress tensor.
 * @param tangent Optional consistent material tangent `dS/dE`.
 */
void NeoHookeElasticity::evaluate_full(const Mat3& C, Mat3& stress, Mat6* tangent) const {
    // Validate the metric before inversion and square-root evaluation.
    const Precision det_c = C.determinant();

    logging::error(det_c > Precision(0),
        "NEO_HOOKE: non-positive right Cauchy-Green determinant");

    // Build the scalar invariants and reusable factors of the constitutive law.
    const Precision J                = std::sqrt(det_c);
    const Precision first_invariant  = C.trace();
    const Precision mean_invariant   = first_invariant / Precision(3);
    const Precision deviatoric_scale = Precision(2) * c10
                                     * std::pow(J, Precision(-2) / Precision(3));
    const Precision volumetric_scale = Precision(2) * J * (J - Precision(1)) / d1;
    const Mat3      C_inv             = C.inverse();

    // Evaluate the isochoric and volumetric PK2 stress contributions.
    stress = deviatoric_scale * (Mat3::Identity() - mean_invariant * C_inv)
           + volumetric_scale * C_inv;

    // A stress-only constitutive query stops here. This is the relevant fast path
    // for residual-only global nonlinear assembly.
    if (tangent == nullptr) {
        return;
    }

    // Differentiate the stress analytically for each independent
    // Green-Lagrange strain component.
    tangent->setZero();
    for (Index column = 0; column < 6; ++column) {
        // Construct one unit Green-Lagrange strain variation. VolumeStrain uses
        // engineering shear components, so tensor shear entries are one half.
        Mat3 dE = Mat3::Zero();

        if (column < 3) {
            dE(column, column) = Precision(1);
        } else if (column == 3) {
            dE(1, 2) = Precision(0.5);
            dE(2, 1) = Precision(0.5);
        } else if (column == 4) {
            dE(0, 2) = Precision(0.5);
            dE(2, 0) = Precision(0.5);
        } else {
            dE(0, 1) = Precision(0.5);
            dE(1, 0) = Precision(0.5);
        }

        // Differentiate C^-1 and the invariants. Because C = I + 2E,
        //
        //     dC = 2 dE
        //     d(C^-1) = -2 C^-1 dE C^-1.
        //
        // The factor two is carried explicitly by the formulas below.
        const Precision inverse_contraction = (C_inv * dE).trace();
        const Precision strain_trace        = dE.trace();
        const Mat3      inverse_derivative  = C_inv * dE * C_inv;

        // Differentiate the isochoric contribution, including the J^(-2/3)
        // scaling, I1/3 term, and inverse metric.
        const Mat3 deviatoric_derivative = deviatoric_scale * (
            -Precision(2) / Precision(3) * inverse_contraction
                * (Mat3::Identity() - mean_invariant * C_inv)
            -Precision(2) / Precision(3) * strain_trace * C_inv
            +Precision(2) * mean_invariant * inverse_derivative
        );

        // Differentiate the volumetric scale and inverse-metric factor.
        const Precision volumetric_derivative = Precision(2) * J
                                              * (Precision(2) * J - Precision(1))
                                              / d1 * inverse_contraction;
        const Mat3 dS = deviatoric_derivative
                      + volumetric_derivative * C_inv
                      - Precision(2) * volumetric_scale * inverse_derivative;

        // Store the physical-stress Voigt derivative corresponding to the current
        // engineering-strain basis direction.
        tangent->col(column) = VolumeStressPK2(dS).voigt();
    }
}

} // namespace fem::material
