/**
 * @file isotropic_elasticity.h
 * @brief Declares homogeneous isotropic linear elasticity.
 *
 * The model supplies constant axial, three-dimensional and shell material
 * tangents from Young's modulus and Poisson's ratio. Green-Lagrange strain and
 * second Piola-Kirchhoff stress use the same constant material tangent.
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
 * The model has no history variables. The state pointers accepted through the
 * common elasticity interface are deliberately unused by every evaluation.
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

    // Total-Lagrangian axial Hooke response S = E E_GL. The optional tangent is
    // the constant derivative dS/dE = E.
    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                old_state,
                  Precision*                      new_state,
                  AxialStressPK2&                 stress,
                  Precision*                      tangent = nullptr) const override;

    // Finite-strain response using the same constant Hooke operator. Input is
    // Green-Lagrange strain and output is second Piola-Kirchhoff stress.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent = nullptr) const override;

    // Green-Lagrange five-component shell response with PK2 output. The material
    // state remains unchanged and the reduced tangent is optional.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

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
