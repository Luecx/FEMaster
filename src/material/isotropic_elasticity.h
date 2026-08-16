/**
 * @file isotropic_elasticity.h
 * @brief Declares homogeneous isotropic linear elasticity.
 *
 * The model supplies constant axial, three-dimensional and shell material
 * tangents from Young's modulus and Poisson's ratio. Linearized calls return
 * Cauchy stress, while Green-Lagrange calls use the identical constant material
 * tangent to return second Piola-Kirchhoff stress.
 *
 * @see Elasticity
 * @see GeneralisedIsotropicElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "elasticity.h"

namespace fem::material {

/**
 * @brief Stateless isotropic Hooke elasticity for axial, solid and shell use.
 *
 * The three-dimensional tangent follows the Lamé form. Shell evaluation uses
 * an in-plane plane-stress block and two transverse shear components. All
 * tangents are expressed in the material basis supplied by the owning section.
 *
 * The model has no history variables. The state pointer accepted through the
 * common elasticity interface is deliberately left unchanged by every
 * evaluation.
 */
struct IsotropicElasticity : Elasticity {
    // Independent elastic constants and derived shear modulus
    Precision youngs;
    Precision poisson;
    Precision shear;

    // Construct the model from Young's modulus and Poisson's ratio. The shear
    // modulus is derived as G = E / (2 (1 + nu)); invalid stability bounds are
    // rejected by the definition in the implementation.
    IsotropicElasticity(Precision youngs_in, Precision poisson_in);

    // Advertise all axial, volume and shell strain measures implemented below.
    // These flags allow the owning section to validate its requested kinematics.
    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Linearized axial Hooke response: sigma = E epsilon and d sigma/d epsilon = E.
    // The supplied state row is valid but unchanged because this law is stateless.
    void evaluate(const AxialStrainLinearized& strain,
                  Precision*                   state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    // Total-Lagrangian axial Hooke response: S = E E_GL. PK2 stress and
    // Green-Lagrange strain are work-conjugate in the reference material axis.
    void evaluate(const AxialStrainGreenLagrange& strain,
                  Precision*                      state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override;

    // Linearized three-dimensional response in material coordinates. The
    // constant isotropic six-by-six tangent maps engineering strain to Cauchy stress.
    void evaluate(const VolumeStrainLinearized& strain,
                  Precision*                    state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    // Finite-strain material response using the same constant Hooke operator.
    // Input is Green-Lagrange strain and output is second Piola-Kirchhoff stress.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  Precision*                       state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    // Linearized five-component shell response with an in-plane plane-stress
    // block and transverse shear moduli. Output stress is Cauchy stress.
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  Precision*                           state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override;

    // Green-Lagrange five-component shell response with PK2 output. The same
    // constant reduced tangent is used and the material state remains unchanged.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  Precision*                              state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    // Build the in-plane plane-stress operator ordered as [11,22,12].
    [[nodiscard]] Mat3 plane_stress_tangent() const;

    // Embed the plane-stress block and transverse shear moduli into the
    // five-component shell material ordering [11,22,12,13,23].
    [[nodiscard]] Mat5 shell_material_tangent() const;

    // Build the full isotropic three-dimensional engineering-Voigt tangent.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
