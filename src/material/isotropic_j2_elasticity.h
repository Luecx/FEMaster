/**
 * @file isotropic_j2_elasticity.h
 * @brief Declares isotropic associative J2 elastoplasticity for axial, solid and shell use.
 *
 * The material supports infinitesimal and finite-strain kinematics. Finite strain
 * uses a multiplicative formulation in which the plastic history is represented
 * by the symmetric plastic metric Cp = Fp^T Fp and the accumulated equivalent
 * plastic strain. The seven stored history scalars are therefore
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp].
 *
 * The public interface follows FEMaster's committed/trial state lifecycle:
 * `old_state` is immutable and every constitutive call starts from that committed
 * history. The return mapping is performed on a private working copy. When
 * `new_state` is non-null the complete converged history is written there; a null
 * `new_state` evaluates exactly the same constitutive response but discards the
 * trial history. This keeps stiffness, prestress and output paths state-neutral
 * without changing the constitutive equations used by nonlinear equilibrium.
 *
 * @see Elasticity
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 02.09.2026
 */

#pragma once

#include "elasticity.h"

#include <algorithm>
#include <array>
#include <vector>

namespace fem::material {

/**
 * @brief Isotropic J2 plasticity with tabulated isotropic hardening.
 *
 * Yield points are supplied as `(yield stress, equivalent plastic strain)` pairs.
 * The first point defines initial yield at zero plastic strain, intermediate
 * points are linearly interpolated and the final yield stress is continued as
 * perfect plasticity beyond the last tabulated strain.
 *
 * Small-strain response uses the standard associative radial-return equations.
 * Finite strain uses the multiplicative plastic metric formulation implemented in
 * `isotropic_j2_elasticity.cpp`. Axial and shell evaluations perform the required
 * transverse/plane-stress reductions around the same three-dimensional update.
 */
struct IsotropicJ2Elasticity : Elasticity {
    struct YieldPoint {
        Precision yield_stress;
        Precision equivalent_plastic_strain;
    };

    Precision youngs;
    Precision poisson;

    IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in);

    [[nodiscard]] Precision shear_modulus() const;
    [[nodiscard]] Precision bulk_modulus() const;

    void add_yield_point(Precision yield_stress, Precision equivalent_plastic_strain);
    [[nodiscard]] const std::vector<YieldPoint>& get_yield_points() const;

    // Supported kinematic/stress-measure pairs.
    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Seven persistent scalars store Cp and accumulated equivalent plastic strain.
    Index state_size() const override;
    void  initialize_state(Precision* state) const override;

    // Constitutive evaluations always integrate from old_state. A non-null
    // new_state receives the converged trial history; nullptr discards it.
    void evaluate(const AxialStrainLinearized& strain,
                  const Precision*             old_state,
                  Precision*                   new_state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                old_state,
                  Precision*                      new_state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              old_state,
                  Precision*                    new_state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     old_state,
                  Precision*                           new_state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override {
        evaluate_from_committed(strain, old_state, new_state, stress, tangent);
    }

private:
    static constexpr Index state_count = 7;

    /**
     * @brief Runs one constitutive update from immutable committed history.
     *
     * The legacy numerical kernels operate in-place on a seven-scalar working
     * row. Keeping that mutation private gives the kernels their natural local
     * representation while enforcing the solver-facing committed/trial contract.
     * Every Newton or line-search candidate therefore starts from the same
     * `old_state`; no previous trial state can feed the next material evaluation.
     */
    template<typename Strain, typename Stress, typename Tangent>
    void evaluate_from_committed(const Strain&    strain,
                                 const Precision* old_state,
                                 Precision*       new_state,
                                 Stress&          stress,
                                 Tangent&         tangent) const {
        std::array<Precision, state_count> working{};
        std::copy_n(old_state, state_count, working.data());

        evaluate(strain, working.data(), stress, tangent);

        if (new_state != nullptr) {
            std::copy_n(working.data(), state_count, new_state);
        }
    }

    // In-place constitutive kernels. They mutate only the private working copy
    // created by evaluate_from_committed(); persistent state ownership stays with
    // the nonlinear state manager.
    void evaluate(const AxialStrainLinearized& strain,
                  Precision*                   state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const;

    void evaluate(const AxialStrainGreenLagrange& strain,
                  Precision*                      state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const;

    void evaluate(const VolumeStrainLinearized& strain,
                  Precision*                    state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const;

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  Precision*                       state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const;

    void evaluate(const ShellMaterialStrainLinearized& strain,
                  Precision*                           state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const;

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  Precision*                              state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const;

    std::vector<YieldPoint> yield_points_;
};

} // namespace fem::material
