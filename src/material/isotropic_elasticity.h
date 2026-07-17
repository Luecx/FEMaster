#pragma once

#include "elasticity.h"

namespace fem::material {

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
                  const Precision*             state_old,
                  Precision*                   state_new,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                state_old,
                  Precision*                      state_new,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override;

    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              state_old,
                  Precision*                    state_new,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 state_old,
                  Precision*                       state_new,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     state_old,
                  Precision*                           state_new,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override;

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        state_old,
                  Precision*                              state_new,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    [[nodiscard]] Mat3 plane_stress_tangent() const;
    [[nodiscard]] Mat5 shell_material_tangent() const;
    [[nodiscard]] Mat6 volume_tangent() const;
};

} // namespace fem::material
