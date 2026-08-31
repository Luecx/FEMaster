/**
 * @file isotropic_j2_elasticity.h
 * @brief Declares isotropic J2 elastoplasticity for axial, solid and shell use.
 *
 * The finite-strain formulation is multiplicative. Plastic history is represented
 * by the symmetric plastic metric Cp = Fp^T Fp together with the accumulated
 * equivalent plastic strain. The total right Cauchy-Green tensor follows from
 * Green-Lagrange strain as C = I + 2 E.
 *
 * Constitutive state is never modified in place. `evaluate()` reads an already
 * established history state and returns the corresponding frozen-history stress
 * and tangent. `integrate()` advances an immutable source state into a distinct
 * target state and is the only operation that performs plastic return mapping.
 *
 * @see Elasticity
 * @see IsotropicElasticity
 */

#pragma once

#include "elasticity.h"

#include <vector>

namespace fem::material {

/**
 * @brief Isotropic associative J2 elastoplasticity with finite-strain support.
 *
 * Hardening points are appended as `(yield stress, equivalent plastic strain)`.
 * The seven scalar history entries are
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp].
 */
struct IsotropicJ2Elasticity : Elasticity {
    struct YieldPoint {
        Precision yield_stress;
        Precision equivalent_plastic_strain;
    };

    Precision youngs;
    Precision poisson;

    [[nodiscard]] Precision shear_modulus() const;
    [[nodiscard]] Precision bulk_modulus() const;

    IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in);

    void add_yield_point(Precision yield_stress, Precision equivalent_plastic_strain);
    [[nodiscard]] const std::vector<YieldPoint>& get_yield_points() const;

    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    Index state_size() const override;
    void  initialize_state(Precision* state) const override;

    // -------------------------------------------------------------------------
    // Read-only response at an already established plastic history state
    // -------------------------------------------------------------------------

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

    // -------------------------------------------------------------------------
    // Plastic integration: immutable source history -> complete target history
    // -------------------------------------------------------------------------

    void integrate(const AxialStrainLinearized& strain,
                   const Precision*             state,
                   Precision*                   target_state,
                   AxialStressCauchy&           stress,
                   Precision&                   tangent) const override;

    void integrate(const AxialStrainGreenLagrange& strain,
                   const Precision*                  state,
                   Precision*                      target_state,
                   AxialStressPK2&                 stress,
                   Precision&                     tangent) const override;

    void integrate(const VolumeStrainLinearized& strain,
                   const Precision*              state,
                   Precision*                    target_state,
                   VolumeStressCauchy&           stress,
                   Mat6&                         tangent) const override;

    void integrate(const VolumeStrainGreenLagrange& strain,
                   const Precision*                   state,
                   Precision*                       target_state,
                   VolumeStressPK2&                 stress,
                   Mat6&                           tangent) const override;

    void integrate(const ShellMaterialStrainLinearized& strain,
                   const Precision*                     state,
                   Precision*                           target_state,
                   ShellMaterialStressCauchy&            stress,
                   Mat5&                               tangent) const override;

    void integrate(const ShellMaterialStrainGreenLagrange& strain,
                   const Precision*                          state,
                   Precision*                              target_state,
                   ShellMaterialStressPK2&                 stress,
                   Mat5&                                  tangent) const override;

private:
    std::vector<YieldPoint> yield_points_;
};

} // namespace fem::material
