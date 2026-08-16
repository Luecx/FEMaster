/**
 * @file constitutive_law.h
 * @brief Declares the solver-facing mechanical constitutive law.
 *
 * A ConstitutiveLaw composes an elastic backbone with an optional plasticity
 * model. It is the only mechanical material interface that element and section
 * code should use for nonlinear constitutive updates.
 *
 * The state contract is transactional:
 *
 *     update(strain, state_old, state_new, stress, tangent)
 *
 * reads committed history from state_old and writes the complete trial history
 * to state_new. The two buffers must not alias for a stateful law. Recovery is
 * state-neutral and never modifies committed or trial history.
 *
 * Existing Elasticity implementations keep their legacy mutable state argument
 * internally for source compatibility. ConstitutiveLaw owns that legacy state
 * slice and prevents it from leaking into the global nonlinear state contract.
 * Current elastic models are stateless, while Plasticity owns the actual history
 * used by J2 and future inelastic laws.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "elasticity.h"
#include "plasticity.h"

#include <memory>

namespace fem::material {

struct ConstitutiveLaw {
    using Ptr = std::shared_ptr<ConstitutiveLaw>;

    bool has_elasticity() const { return elasticity_ != nullptr; }
    bool has_plasticity() const { return plasticity_ != nullptr; }

    const ElasticityPtr& elasticity() const { return elasticity_; }
    const PlasticityPtr& plasticity() const { return plasticity_; }

    void set_elasticity(ElasticityPtr elasticity);
    void set_plasticity(PlasticityPtr plasticity);

    Index elastic_state_size() const;
    Index plastic_state_size() const;
    Index state_size() const;
    void initialize_state(Precision* state) const;

    bool supports_axial_linearized() const;
    bool supports_axial_green_lagrange() const;
    bool supports_volume_linearized() const;
    bool supports_volume_green_lagrange() const;
    bool supports_beam_resultants() const;
    bool supports_shell_integration_linearized() const;
    bool supports_shell_integration_green_lagrange() const;

    void update(const AxialStrainLinearized& strain,
                const Precision* state_old,
                Precision* state_new,
                AxialStressCauchy& stress,
                Precision& tangent) const;
    void update(const AxialStrainGreenLagrange& strain,
                const Precision* state_old,
                Precision* state_new,
                AxialStressPK2& stress,
                Precision& tangent) const;
    void update(const VolumeStrainLinearized& strain,
                const Precision* state_old,
                Precision* state_new,
                VolumeStressCauchy& stress,
                Mat6& tangent) const;
    void update(const VolumeStrainGreenLagrange& strain,
                const Precision* state_old,
                Precision* state_new,
                VolumeStressPK2& stress,
                Mat6& tangent) const;
    void update(const BeamGeneralizedStrain& strain,
                const Precision* state_old,
                Precision* state_new,
                BeamStressResultants& resultants,
                Mat6& tangent) const;
    void update(const ShellMaterialStrainLinearized& strain,
                const Precision* state_old,
                Precision* state_new,
                ShellMaterialStressCauchy& stress,
                Mat5& tangent) const;
    void update(const ShellMaterialStrainGreenLagrange& strain,
                const Precision* state_old,
                Precision* state_new,
                ShellMaterialStressPK2& stress,
                Mat5& tangent) const;

    void recover(const AxialStrainLinearized& strain,
                 const Precision* state,
                 AxialStressCauchy& stress) const;
    void recover(const AxialStrainGreenLagrange& strain,
                 const Precision* state,
                 AxialStressPK2& stress) const;
    void recover(const VolumeStrainLinearized& strain,
                 const Precision* state,
                 VolumeStressCauchy& stress) const;
    void recover(const VolumeStrainGreenLagrange& strain,
                 const Precision* state,
                 VolumeStressPK2& stress) const;
    void recover(const ShellMaterialStrainLinearized& strain,
                 const Precision* state,
                 ShellMaterialStressCauchy& stress) const;
    void recover(const ShellMaterialStrainGreenLagrange& strain,
                 const Precision* state,
                 ShellMaterialStressPK2& stress) const;

private:
    void validate() const;
    Precision* prepare_elastic_trial_state(const Precision* state_old,
                                           Precision* state_new) const;
    const Precision* plastic_state_old(const Precision* state_old) const;
    Precision* plastic_state_new(Precision* state_new) const;

    ElasticityPtr elasticity_ = nullptr;
    PlasticityPtr plasticity_ = nullptr;
};

using ConstitutiveLawPtr = ConstitutiveLaw::Ptr;

} // namespace fem::material
