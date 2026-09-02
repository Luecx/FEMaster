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
 * contract. `old_state` is read-only and `new_state` receives the converged
 * candidate history when non-null. A null tangent requests the same physical
 * constitutive update without material-tangent output. The full-volume J2 path
 * can therefore skip its expensive algorithmic tangent during residual-only
 * assembly; axial and shell reductions still require local derivatives to solve
 * their traction-free constitutive constraints.
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
 * values are linearly interpolated. Beyond the final point the final yield
 * stress is continued, so the material becomes perfectly plastic after the
 * tabulated range.
 *
 * The small-strain model uses the classical associative radial return. The
 * finite-strain model uses the multiplicative split
 *
 *     F = Fe Fp
 *
 * together with the exponential plastic update
 *
 *     Fp_(n+1) = exp(A) Fp_n.
 *
 * Its algorithmic tangent is the consistent analytic derivative of the same
 * discrete return map used for the stress update. Persistent state contains only
 * `Cp` and accumulated equivalent plastic strain; local return-map unknowns are
 * reconstructed for every candidate from the committed state.
 */
struct IsotropicJ2Elasticity : Elasticity {
    using Elasticity::evaluate;

    /**
     * @brief One point of the piecewise-linear isotropic hardening law.
     *
     * `equivalent_plastic_strain` is the accumulated scalar J2 history variable
     * and `yield_stress` is the corresponding current radius of the yield
     * surface. Points are stored in strictly increasing plastic-strain order.
     */
    struct YieldPoint {
        Precision yield_stress;
        Precision equivalent_plastic_strain;
    };

    // Independent elastic parameters
    Precision youngs;
    Precision poisson;

    // Construction and elastic constants
    IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in);

    [[nodiscard]] Precision shear_modulus() const;
    [[nodiscard]] Precision bulk_modulus() const;

    // Piecewise-linear isotropic hardening definition. The first point must be
    // located at zero equivalent plastic strain; subsequent points must increase
    // monotonically in both plastic strain and yield stress.
    void add_yield_point(Precision yield_stress, Precision equivalent_plastic_strain);
    [[nodiscard]] const std::vector<YieldPoint>& get_yield_points() const;

    // Supported strain/stress-measure pairs
    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Persistent material-point history. Seven scalars store the six independent
    // components of Cp = Fp^T Fp followed by accumulated equivalent plastic strain.
    Index state_size() const override;
    void  initialize_state(Precision* state) const override;

    // Axial reductions enforce zero transverse stress through the corresponding
    // three-dimensional constitutive law. The external scalar tangent is optional,
    // although the finite/local reduction may require derivatives internally.
    void evaluate(const AxialStrainLinearized& strain,
                  const Precision*             old_state,
                  Precision*                   new_state,
                  AxialStressCauchy&           stress,
                  Precision*                   tangent = nullptr) const override;

    void evaluate(const AxialStrainGreenLagrange& strain,
                  const Precision*                old_state,
                  Precision*                      new_state,
                  AxialStressPK2&                 stress,
                  Precision*                      tangent = nullptr) const override;

    // Three-dimensional material response. A null tangent preserves the complete
    // stress/state update while omitting construction of the algorithmic matrix.
    void evaluate(const VolumeStrainLinearized& strain,
                  const Precision*              old_state,
                  Precision*                    new_state,
                  VolumeStressCauchy&           stress,
                  Mat6*                         tangent = nullptr) const override;

    void evaluate(const VolumeStrainGreenLagrange& strain,
                  const Precision*                 old_state,
                  Precision*                       new_state,
                  VolumeStressPK2&                 stress,
                  Mat6*                            tangent = nullptr) const override;

    // Integrated-shell material response under the local plane-stress condition
    // S33 = 0. Tangent output is optional; the thickness solve itself still uses
    // the three-dimensional consistent derivative for its Newton update.
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  const Precision*                     old_state,
                  Precision*                           new_state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5*                                tangent = nullptr) const override;

    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;

private:
    // Persistent state layout size
    static constexpr Index state_count = 7;

    // Piecewise-linear isotropic hardening table
    std::vector<YieldPoint> yield_points_;
};

} // namespace fem::material
