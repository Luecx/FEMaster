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
 * The model has no history variables. The source state accepted through the
 * common elasticity interface is read-only and ignored by every evaluation.
 */
struct IsotropicElasticity : Elasticity {
    Precision youngs;
    Precision poisson;
    Precision shear;

    IsotropicElasticity(Precision youngs_in, Precision poisson_in);

    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    void evaluate(const AxialStrainLinearized& strain,
                  const Precision*             state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                  state,
                  AxialStressPK2&                 stress,
                  Precision&                     tangent) const override;

    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                   state,
                  VolumeStressPK2&                 stress,
                  Mat6&                           tangent) const override;

    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                               tangent) const override;

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                          state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                  tangent) const override;

private:
    [[nodiscard]] Mat3 plane_stress_tangent() const;
    [[nodiscard]] Mat5 shell_material_tangent() const;
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
