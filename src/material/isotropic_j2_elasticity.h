/**
 * @file isotropic_j2_elasticity.h
 * @brief Declares isotropic associative J2 elastoplasticity.
 *
 * The material supports both infinitesimal and finite-strain kinematics. The
 * persistent constitutive history contains only physical state variables:
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp]
 *
 * where Cp = Fp^T Fp is the symmetric plastic metric and eqp is the accumulated
 * equivalent plastic strain. Local Newton unknowns used by the finite-strain
 * return map are deliberately not persisted between global increments.
 *
 * Every public constitutive evaluation follows FEMaster's committed/trial state
 * contract. `old_state` is read-only, `new_state` receives the converged candidate
 * history when non-null, and the same physical update is evaluated when
 * `new_state == nullptr` without modifying persistent history.
 *
 * @see Elasticity
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 02.09.2026
 */

#pragma once

#include "elasticity.h"

#include <vector>

namespace fem::material {

/**
 * @brief Isotropic J2 plasticity with tabulated isotropic hardening.
 *
 * Yield points are supplied as `(yield stress, equivalent plastic strain)` pairs.
 * The first point defines initial yield at zero plastic strain. Intermediate
 * values are linearly interpolated. Beyond the final point the final yield stress
 * is continued, i.e. the material is perfectly plastic after the tabulated range.
 *
 * The small-strain model uses the classical associative radial return. The
 * finite-strain model uses a multiplicative split F = Fe Fp with an exponential
 * plastic update and an analytically consistent algorithmic tangent.
 */
struct IsotropicJ2Elasticity : Elasticity {
    using Elasticity::evaluate;

    /**
     * @brief One point of the piecewise-linear isotropic hardening curve.
     */
    struct YieldPoint {
        Precision yield_stress;
        Precision equivalent_plastic_strain;
    };

    Precision youngs;
    Precision poisson;

    /**
     * @brief Constructs isotropic J2 plasticity from Young's modulus and Poisson ratio.
     * @param youngs_in Young's modulus E.
     * @param poisson_in Poisson ratio nu.
     */
    IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in);

    /**
     * @brief Returns the elastic shear modulus G = E / (2 (1 + nu)).
     */
    [[nodiscard]] Precision shear_modulus() const;

    /**
     * @brief Returns the elastic bulk modulus K = E / (3 (1 - 2 nu)).
     */
    [[nodiscard]] Precision bulk_modulus() const;

    /**
     * @brief Appends one point to the isotropic hardening curve.
     * @param yield_stress Current uniaxial yield stress.
     * @param equivalent_plastic_strain Accumulated equivalent plastic strain.
     */
    void add_yield_point(Precision yield_stress, Precision equivalent_plastic_strain);

    /**
     * @brief Returns the complete tabulated isotropic hardening curve.
     */
    [[nodiscard]] const std::vector<YieldPoint>& get_yield_points() const;

    /** @brief Reports support for infinitesimal axial strain with Cauchy stress. */
    bool supports_axial_linearized() const override;

    /** @brief Reports support for Green-Lagrange axial strain with PK2 stress. */
    bool supports_axial_green_lagrange() const override;

    /** @brief Reports support for infinitesimal three-dimensional strain. */
    bool supports_volume_linearized() const override;

    /** @brief Reports support for Green-Lagrange three-dimensional strain. */
    bool supports_volume_green_lagrange() const override;

    /** @brief Reports support for infinitesimal integrated-shell material response. */
    bool supports_shell_integration_linearized() const override;

    /** @brief Reports support for finite-strain integrated-shell material response. */
    bool supports_shell_integration_green_lagrange() const override;

    /**
     * @brief Returns the number of persistent constitutive state scalars.
     *
     * The seven values are the six independent components of Cp and eqp.
     */
    Index state_size() const override;

    /**
     * @brief Initializes the constitutive state to Fp = I and eqp = 0.
     * @param state Writable state row containing at least state_size() entries.
     */
    void initialize_state(Precision* state) const override;

    /**
     * @brief Evaluates the infinitesimal uniaxial J2 response.
     */
    void evaluate(const AxialStrainLinearized& strain,
                  const Precision*             old_state,
                  Precision*                   new_state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    /**
     * @brief Evaluates the finite-strain uniaxial J2 response.
     */
    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                old_state,
                  Precision*                      new_state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override;

    /**
     * @brief Evaluates the infinitesimal three-dimensional J2 response.
     */
    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              old_state,
                  Precision*                    new_state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    /**
     * @brief Evaluates the finite-strain three-dimensional J2 response.
     */
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    /**
     * @brief Evaluates the finite-strain three-dimensional response with optional tangent.
     *
     * A null tangent requests the identical return map and trial state but skips
     * construction of the algorithmic material tangent. This is used by residual-
     * only global line-search evaluations.
     */
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent) const override;

    /**
     * @brief Evaluates the infinitesimal integrated-shell J2 response under S33 = 0.
     */
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     old_state,
                  Precision*                           new_state,
                  ShellMaterialStressCauchy&           stress,
                  Mat5&                                tangent) const override;

    /**
     * @brief Evaluates the finite-strain integrated-shell J2 response under S33 = 0.
     */
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    static constexpr Index state_count = 7;

    std::vector<YieldPoint> yield_points_;
};

} // namespace fem::material
