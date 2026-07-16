#pragma once

#include "elasticity.h"

namespace fem::material {

struct OrthotropicElasticity : Elasticity {
    Precision Ex;
    Precision Ey;
    Precision Ez;
    Precision Gyz;
    Precision Gzx;
    Precision Gxy;
    Precision vyz;
    Precision vzx;
    Precision vxy;

    OrthotropicElasticity(Precision ex,
                          Precision ey,
                          Precision ez,
                          Precision gyz,
                          Precision gzx,
                          Precision gxy,
                          Precision vyz,
                          Precision vzx,
                          Precision vxy);

    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_resultants() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

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

    void evaluate(const ShellGeneralizedStrain& strain,
                  Precision                     thickness,
                  const Precision*              state_old,
                  Precision*                    state_new,
                  ShellStressResultants&        resultants,
                  Mat8&                         tangent) const override;

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
    [[nodiscard]] Mat8 shell_resultant_tangent(Precision thickness) const;
};

} // namespace fem::material
