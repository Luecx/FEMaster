/**
 * @file orthotropic_elasticity.h
 * @brief Declares homogeneous orthotropic linear elasticity.
 *
 * The material is defined by the nine conventional engineering constants
 * `E1`, `E2`, `E3`, `nu12`, `nu13`, `nu23`, `G12`, `G13` and `G23` in its local
 * material basis. This parameterization matches the standard orthotropic
 * engineering-constants representation used by Abaqus and avoids storing
 * reciprocal Poisson ratios as independent material data.
 *
 * Solid and shell sections supply the coordinate transformation between the
 * local material basis and the element basis. The three-dimensional material
 * tangent is recovered from the symmetric engineering compliance, while shell
 * response uses the corresponding orthotropic plane-stress reduction.
 *
 * @see Elasticity
 * @see SolidSection
 * @see IntegratedShellSection
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "elasticity.h"

namespace fem::material {

/**
 * @brief Stateless orthotropic Hooke elasticity for solid and shell response.
 *
 * The stored material parameters are the three directional Young's moduli,
 * three major Poisson ratios and three engineering shear moduli in principal
 * material directions. Reciprocal Poisson ratios follow from symmetry of the
 * compliance matrix and are therefore derived rather than stored.
 *
 * The volume tangent is obtained by inverting the symmetric engineering
 * compliance. Shell response uses an in-plane orthotropic plane-stress tangent
 * together with the independent `G13` and `G23` transverse shear moduli.
 * Linearized evaluation returns Cauchy stress; Green-Lagrange evaluation returns
 * second Piola-Kirchhoff stress.
 */
struct OrthotropicElasticity : Elasticity {
    // Young's moduli along the principal material directions
    Precision E1;
    Precision E2;
    Precision E3;

    // Major Poisson ratios nu_ij: transverse strain in j for loading in i
    Precision nu12;
    Precision nu13;
    Precision nu23;

    // Engineering shear moduli in the principal material planes
    Precision G12;
    Precision G13;
    Precision G23;

    // Construct the material directly from conventional orthotropic engineering
    // constants. Reciprocal Poisson ratios are derived from material symmetry.
    OrthotropicElasticity(Precision E1,
                          Precision E2,
                          Precision E3,
                          Precision nu12,
                          Precision nu13,
                          Precision nu23,
                          Precision G12,
                          Precision G13,
                          Precision G23);

    // Advertise three-dimensional and shell response for both infinitesimal and
    // Green-Lagrange strain measures. Axial and beam reductions are unsupported.
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Linearized three-dimensional orthotropic response in material axes.
    // tangent maps engineering strain to Cauchy stress and state is unchanged.
    void evaluate(const VolumeStrainLinearized& strain,
                  Precision*                    state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    // Total-Lagrangian orthotropic response using the same constant material
    // operator, interpreted as dS/dE for PK2 stress and Green-Lagrange strain.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  Precision*                       state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    // Linearized shell response in the material basis. In-plane terms use the
    // orthotropic plane-stress reduction; 13 and 23 shear use G13 and G23.
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  Precision*                           state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override;

    // Finite-strain shell response returning PK2 components work-conjugate to
    // the five supplied Green-Lagrange strain components. State remains unchanged.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  Precision*                              state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    // Build the in-plane orthotropic plane-stress tangent ordered [11,22,12].
    [[nodiscard]] Mat3 plane_stress_tangent() const;

    // Embed the in-plane block and directional transverse shear moduli into the
    // five-component shell material ordering.
    [[nodiscard]] Mat5 shell_material_tangent() const;

    // Invert the symmetric engineering compliance to obtain the full
    // three-dimensional tangent in local material coordinates.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
