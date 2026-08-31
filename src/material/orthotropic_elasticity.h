/**
 * @file orthotropic_elasticity.h
 * @brief Declares homogeneous orthotropic linear elasticity.
 *
 * The material is defined by the nine conventional engineering constants
 * `E1`, `E2`, `E3`, `nu12`, `nu13`, `nu23`, `G12`, `G13` and `G23` in its local
 * material basis.
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
 */
struct OrthotropicElasticity : Elasticity {
    Precision E1;
    Precision E2;
    Precision E3;

    Precision nu12;
    Precision nu13;
    Precision nu23;

    Precision G12;
    Precision G13;
    Precision G23;

    OrthotropicElasticity(Precision E1,
                          Precision E2,
                          Precision E3,
                          Precision nu12,
                          Precision nu13,
                          Precision nu23,
                          Precision G12,
                          Precision G13,
                          Precision G23);

    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

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
