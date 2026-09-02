/**
 * @file neo_hooke_elasticity.h
 * @brief Declares compressible isotropic Neo-Hookean elasticity.
 *
 * The finite-strain response follows the decoupled Abaqus-type potential
 *
 *     W = C10 (I1_bar - 3) + (J - 1)^2 / D1.
 *
 * Solid evaluation returns second Piola-Kirchhoff stress and, when requested,
 * its consistent material tangent. Axial and shell reductions enforce
 * traction-free lateral or thickness response through local Newton iterations.
 * Linearized overloads use the infinitesimal tangent derived from the same
 * material parameters.
 *
 * @see Elasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "elasticity.h"

namespace fem::material {

/**
 * @brief Stateless compressible isotropic Neo-Hookean material.
 *
 * The three-dimensional formulation is expressed through the right Cauchy-
 * Green tensor in the reference material basis. Its constitutive output is
 * second Piola-Kirchhoff stress and, optionally, the consistent derivative with
 * respect to Green-Lagrange strain in engineering-Voigt ordering.
 *
 * The axial reduction solves for two equal transverse stretches such that the
 * lateral PK2 stresses vanish. The shell reduction solves the missing thickness
 * component of the right Cauchy-Green tensor such that `S33 = 0`. These local
 * solves require the full material derivative internally even for a stress-only
 * query, whereas the full three-dimensional response can skip tangent assembly
 * entirely when the caller passes a null tangent pointer.
 *
 * The local iterations do not represent material history, so the common
 * material-point state remains unchanged.
 */
struct NeoHookeElasticity : Elasticity {
    // Input potential parameters and derived infinitesimal moduli
    Precision c10;
    Precision d1;
    Precision mu{};
    Precision bulk{};
    Precision lame_lambda{};

    // Construct the model from the isochoric coefficient C10 and volumetric
    // penalty parameter D1. The infinitesimal shear, bulk and Lamé moduli are
    // derived once and reused by the linearized response and initial guesses.
    NeoHookeElasticity(Precision c10_in, Precision d1_in);

    // Advertise axial, full-volume and shell reductions for both infinitesimal
    // and Green-Lagrange kinematics.
    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Infinitesimal axial response using the Young's modulus implied by the
    // Neo-Hookean parameters. Stress is Cauchy stress and state is unchanged.
    void evaluate(const AxialStrainLinearized& strain,
                  const Precision*             old_state,
                  Precision*                   new_state,
                  AxialStressCauchy&           stress,
                  Precision*                   tangent = nullptr) const override;

    // Finite axial response reduced from the full three-dimensional potential.
    // The local traction-free solve requires a full tangent internally; the
    // condensed axial tangent is returned only when requested.
    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                old_state,
                  Precision*                      new_state,
                  AxialStressPK2&                 stress,
                  Precision*                      tangent = nullptr) const override;

    // Infinitesimal three-dimensional response using the linearization of the
    // finite-strain potential at C = I. Output is Cauchy stress in material axes.
    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              old_state,
                  Precision*                    new_state,
                  VolumeStressCauchy&           stress,
                  Mat6*                         tangent = nullptr) const override;

    // Full finite-strain response. Green-Lagrange strain is mapped to C = I+2E;
    // the expensive analytic tangent is assembled only when requested.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent = nullptr) const override;

    // Infinitesimal shell response from the linearized plane-stress and
    // transverse-shear tangent implied by the three-dimensional parameters.
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     old_state,
                  Precision*                           new_state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5*                                tangent = nullptr) const override;

    // Finite-strain shell response. The thickness metric is solved so S33 = 0;
    // the condensed tangent is copied to the caller only when requested.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

private:
    // Evaluate full PK2 stress from a positive-definite right Cauchy-Green
    // tensor. The analytic dS/dE matrix is assembled only for a non-null tangent.
    void evaluate_full(const Mat3& right_cauchy_green,
                       Mat3&       stress,
                       Mat6*       tangent = nullptr) const;

    // Linearize the full potential at the undeformed configuration in
    // three-dimensional engineering-Voigt ordering.
    [[nodiscard]] Mat6 linear_tangent() const;

    // Build the infinitesimal five-component shell reduction used by linear shells.
    [[nodiscard]] Mat5 linear_shell_tangent() const;
};

} // namespace fem::material
