/**
 * @file isotropic_j2_elasticity.h
 * @brief Declares isotropic J2 elastoplasticity for axial, solid and shell use.
 *
 * The finite-strain formulation is multiplicative. Plastic history is represented
 * by the symmetric plastic metric Cp = Fp^T Fp together with the accumulated
 * equivalent plastic strain. The total right Cauchy-Green tensor is reconstructed
 * from Green-Lagrange strain as C = I + 2 E, so the isotropic constitutive law
 * does not require the deformation gradient in the public material interface.
 *
 * Finite elasticity uses the same Saint-Venant-Kirchhoff elastic law as
 * IsotropicElasticity, applied to the elastic Green-Lagrange strain in the
 * intermediate configuration. Yielding is evaluated from the deviatoric Mandel
 * stress, the flow rule is associative, and the plastic update uses an exponential
 * map. Isotropic hardening is defined by yield points that are appended to the
 * material after construction. This allows parsers to build the hardening curve
 * incrementally without knowing the number of points in advance.
 *
 * Shell calls enforce the Total-Lagrangian plane-stress convention S33 = 0 by a
 * local constitutive Newton solve and static condensation. Axial calls analogously
 * eliminate the two transverse strains by S22 = S33 = 0.
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
 * Hardening is added incrementally as ordered pairs of uniaxial yield stress and
 * accumulated equivalent plastic strain, matching the common tabular plasticity
 * convention:
 *
 *     add_yield_point(yield_stress, equivalent_plastic_strain)
 *
 * The first point must have zero equivalent plastic strain. Every following point
 * must increase the equivalent plastic strain and may keep or increase the yield
 * stress. Between points the yield stress is interpolated linearly. Beyond the
 * final point it remains constant, i.e. the response becomes perfectly plastic
 * until another point is appended.
 *
 * The seven scalar history entries are
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp]
 *
 * where Cp starts as the identity tensor and eqp is the accumulated equivalent
 * plastic strain. The same history representation is used by linearized and
 * finite-strain evaluations.
 */
struct IsotropicJ2Elasticity : Elasticity {
    struct YieldPoint {
        Precision yield_stress;
        Precision equivalent_plastic_strain;
    };

    // Elastic constants. The hardening curve is populated afterwards through
    // add_yield_point(), so the number of plastic data points need not be known
    // when the material object is constructed.
    Precision youngs;
    Precision poisson;

    // Derived isotropic moduli. They are functions rather than independent
    // material parameters so they cannot drift away from E and nu.
    [[nodiscard]] Precision shear_modulus() const;
    [[nodiscard]] Precision bulk_modulus() const;

    // Construct the elastic part. At least one yield point must be added before
    // the material is evaluated.
    IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in);

    // Append one isotropic-hardening point in the conventional order
    // (yield stress, accumulated equivalent plastic strain). The first point must
    // be at zero plastic strain; subsequent points must be strictly increasing in
    // plastic strain and non-decreasing in yield stress.
    void add_yield_point(Precision yield_stress, Precision equivalent_plastic_strain);

    // Read-only access is useful for input/output code without bypassing the
    // validation performed by add_yield_point().
    [[nodiscard]] const std::vector<YieldPoint>& get_yield_points() const;

    // The model supports all existing axial, volume and integrated-shell paths.
    bool supports_axial_linearized() const override;
    bool supports_axial_green_lagrange() const override;
    bool supports_volume_linearized() const override;
    bool supports_volume_green_lagrange() const override;
    bool supports_shell_integration_linearized() const override;
    bool supports_shell_integration_green_lagrange() const override;

    // Cp has six independent symmetric components plus equivalent plastic strain.
    Index state_size() const override;
    void  initialize_state(Precision* state) const override;

    // Small-strain uniaxial-stress J2 response.
    void evaluate(const AxialStrainLinearized& strain,
                  Precision*                   state,
                  AxialStressCauchy&           stress,
                  Precision&                   tangent) const override;

    // Finite-strain uniaxial-stress J2 response. The two transverse
    // Green-Lagrange strains are eliminated by S22 = S33 = 0.
    void evaluate(const AxialStrainGreenLagrange& strain,
                  Precision*                      state,
                  AxialStressPK2&                 stress,
                  Precision&                      tangent) const override;

    // Small-strain three-dimensional J2 response.
    void evaluate(const VolumeStrainLinearized& strain,
                  Precision*                    state,
                  VolumeStressCauchy&           stress,
                  Mat6&                         tangent) const override;

    // Multiplicative finite-strain J2 response. Input E reconstructs C = I + 2E;
    // output is reference PK2 stress and dS/dE for Total-Lagrangian assembly.
    void evaluate(const VolumeStrainGreenLagrange& strain,
                  Precision*                       state,
                  VolumeStressPK2&                 stress,
                  Mat6&                            tangent) const override;

    // Small-strain shell response with sigma33 = 0.
    void evaluate(const ShellMaterialStrainLinearized& strain,
                  Precision*                           state,
                  ShellMaterialStressCauchy&            stress,
                  Mat5&                                 tangent) const override;

    // Finite-strain shell response with the existing TL plane-stress condition
    // S33 = 0 and a consistently condensed five-component material tangent.
    void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                  Precision*                              state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5&                                   tangent) const override;

private:
    std::vector<YieldPoint> yield_points_;
};

} // namespace fem::material
