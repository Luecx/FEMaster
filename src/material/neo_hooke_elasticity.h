#pragma once

#include "elasticity.h"

namespace fem::material {

/**
 * Compressible isotropic Neo-Hookean elasticity using the Abaqus potential.
 *
 *   W = C10 (I1_bar - 3) + (J - 1)^2 / D1
 *
 * Truss:
 *   Three-dimensional model reduced under
 *   traction-free lateral surfaces.
 */
struct NeoHookeElasticity : Elasticity {
    Precision c10;
    Precision d1;
    Precision mu{};
    Precision bulk{};
    Precision lame_lambda{};

    NeoHookeElasticity(Precision c10_in, Precision d1_in);

    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    void evaluate(const AxialStrainLinearized& strain,
                  Precision*                   state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    void evaluate(const AxialStrainGreenLagrange& strain,
                  Precision*                      state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override;

    void evaluate(const VolumeStrainLinearized& strain,
                  Precision*                    state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  Precision*                       state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    void evaluate(const ShellMaterialStrainLinearized& strain,
                  Precision*                           state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override;

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  Precision*                              state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    void               evaluate_full(const Mat3& right_cauchy_green, Mat3& stress, Mat6& tangent) const;
    [[nodiscard]] Mat6 linear_tangent() const;
    [[nodiscard]] Mat5 linear_shell_tangent() const;
};

} // namespace fem::material
