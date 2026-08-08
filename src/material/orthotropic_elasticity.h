/**
 * @file orthotropic_elasticity.h
 * @brief Declares homogeneous orthotropic linear elasticity.
 *
 * The material is defined by three Young's moduli, three engineering shear
 * moduli and three major Poisson ratios in its local material basis. Solid and
 * shell sections supply the coordinate transformation between that basis and
 * the element basis.
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
 * The three-dimensional tangent is obtained by inverting the symmetric
 * engineering compliance matrix with reciprocal Poisson ratios. Shell response
 * uses the in-plane orthotropic plane-stress tangent together with the two
 * transverse shear moduli. Linearized evaluation returns Cauchy stress;
 * Green-Lagrange evaluation returns second Piola-Kirchhoff stress.
 */
struct OrthotropicElasticity : Elasticity {
    // Young's moduli along the local material axes
    Precision Ex;
    Precision Ey;
    Precision Ez;

    // Engineering shear moduli in the local material planes
    Precision Gyz;
    Precision Gzx;
    Precision Gxy;

    // Major Poisson ratios
    Precision vyz;
    Precision vzx;
    Precision vxy;

    // Construct the local engineering compliance from directional Young's and
    // shear moduli plus the three major Poisson ratios. Reciprocal minor ratios
    // are derived from symmetry when the volume tangent is built.
    OrthotropicElasticity(Precision ex,
                          Precision ey,
                          Precision ez,
                          Precision gyz,
                          Precision gzx,
                          Precision gxy,
                          Precision vyz,
                          Precision vzx,
                          Precision vxy);

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
    // orthotropic plane-stress reduction; XZ and YZ use Gzx and Gyz respectively.
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

    // Invert the reciprocal, symmetric engineering compliance to obtain the
    // full three-dimensional tangent in local material coordinates.
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
