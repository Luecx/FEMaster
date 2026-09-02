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
#include "../core/logging.h"
#include "strain/axial_strain_green_lagrange.h"
#include "stress/axial_stress_pk2.h"

#include <algorithm>
#include <array>
#include <cmath>
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
        log_axial_finite_debug(strain, old_state, new_state, stress, tangent);
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

    /**
     * Emits focused diagnostics for the finite-strain axial J2 path.
     *
     * The nonlinear solver normally suppresses element/material logging during
     * assembly. This debug branch temporarily re-enables logging only while the
     * diagnostic block is emitted, restoring the previous logging state
     * immediately afterwards. Besides the actual PK2 stress and algorithmic
     * tangent, the routine independently finite-differences both dS/dE and the
     * nominal axial response P=lambda*S from the same committed state. The latter
     * is the exact scalar quantity entering the Total-Lagrangian truss residual,
     * so its derivative must match S + lambda^2*dS/dE.
     */
    void log_axial_finite_debug(const AxialStrainGreenLagrange& strain,
                                const Precision*                old_state,
                                const Precision*                new_state,
                                const AxialStressPK2&           stress,
                                Precision                       tangent) const {
        const Precision E = strain.value();
        const Precision lambda_sq = Precision(1) + Precision(2) * E;
        if (!(lambda_sq > Precision(0))) {
            return;
        }

        const Precision lambda = std::sqrt(lambda_sq);
        const Precision h_E = Precision(1e-7) * std::max(Precision(1), std::abs(E));
        const Precision h_lambda = Precision(1e-7) * std::max(Precision(1), std::abs(lambda));

        auto evaluate_stress_from_committed = [&](Precision E_probe) {
            std::array<Precision, state_count> working{};
            std::copy_n(old_state, state_count, working.data());

            AxialStressPK2 probe_stress;
            Precision      probe_tangent = Precision(0);
            evaluate(
                AxialStrainGreenLagrange(E_probe),
                working.data(),
                probe_stress,
                probe_tangent
            );
            return probe_stress.value();
        };

        const Precision S_plus_E  = evaluate_stress_from_committed(E + h_E);
        const Precision S_minus_E = evaluate_stress_from_committed(E - h_E);
        const Precision tangent_fd = (S_plus_E - S_minus_E) / (Precision(2) * h_E);

        const Precision lambda_plus  = lambda + h_lambda;
        const Precision lambda_minus = lambda - h_lambda;
        const Precision E_plus_lambda  = Precision(0.5) * (lambda_plus  * lambda_plus  - Precision(1));
        const Precision E_minus_lambda = Precision(0.5) * (lambda_minus * lambda_minus - Precision(1));
        const Precision S_plus_lambda  = evaluate_stress_from_committed(E_plus_lambda);
        const Precision S_minus_lambda = evaluate_stress_from_committed(E_minus_lambda);

        const Precision nominal_stress = lambda * stress.value();
        const Precision nominal_tangent = stress.value() + lambda * lambda * tangent;
        const Precision nominal_tangent_fd = (
            lambda_plus * S_plus_lambda - lambda_minus * S_minus_lambda
        ) / (Precision(2) * h_lambda);

        const Precision tangent_scale = std::max(Precision(1), std::abs(tangent_fd));
        const Precision nominal_scale = std::max(Precision(1), std::abs(nominal_tangent_fd));
        const Precision tangent_rel_error = std::abs(tangent - tangent_fd) / tangent_scale;
        const Precision nominal_rel_error = std::abs(nominal_tangent - nominal_tangent_fd) / nominal_scale;

        const bool logging_was_enabled = logging::is_enabled();
        if (!logging_was_enabled) {
            logging::enable();
        }

        logging::info(true, "J2DBG ----------------------------------------------------------------");
        logging::info(true, "J2DBG axial finite: E_GL=", E,
                            ", stretch=", lambda,
                            ", S_PK2=", stress.value(),
                            ", P=lambda*S=", nominal_stress);
        logging::info(true, "J2DBG tangent: C_alg=", tangent,
                            ", dS/dE_FD=", tangent_fd,
                            ", rel_err=", tangent_rel_error);
        logging::info(true, "J2DBG truss scalar tangent: S+lambda^2*C=", nominal_tangent,
                            ", d(lambda*S)/dlambda_FD=", nominal_tangent_fd,
                            ", rel_err=", nominal_rel_error);
        logging::info(true, "J2DBG committed state: Cp11=", old_state[0],
                            ", Cp22=", old_state[1],
                            ", Cp33=", old_state[2],
                            ", Cp23=", old_state[3],
                            ", Cp13=", old_state[4],
                            ", Cp12=", old_state[5],
                            ", eqp=", old_state[6]);

        if (new_state != nullptr) {
            logging::info(true, "J2DBG trial state:     Cp11=", new_state[0],
                                ", Cp22=", new_state[1],
                                ", Cp33=", new_state[2],
                                ", Cp23=", new_state[3],
                                ", Cp13=", new_state[4],
                                ", Cp12=", new_state[5],
                                ", eqp=", new_state[6]);
            logging::info(true, "J2DBG state delta:     dCp11=", new_state[0] - old_state[0],
                                ", dCp22=", new_state[1] - old_state[1],
                                ", dCp33=", new_state[2] - old_state[2],
                                ", deqp=", new_state[6] - old_state[6]);
        } else {
            logging::info(true, "J2DBG trial state: discarded (new_state=nullptr)");
        }

        logging::info(true, "J2DBG FD steps: h_E=", h_E, ", h_lambda=", h_lambda);

        if (!logging_was_enabled) {
            logging::disable();
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
