/**
 * @file elasticity.h
 * @brief Declares the common interface for elastic constitutive models.
 *
 * `Elasticity` separates element and section kinematics from constitutive
 * evaluation. Callers select the overload matching their strain measure and
 * dimensional reduction; implementations return the work-conjugate stress and
 * consistent tangent in the same material basis.
 *
 * Every evaluation receives separate input and output state rows for the current
 * material point. The input row is always valid and remains unchanged, while a
 * history-dependent implementation writes its updated history to the output
 * row. Ownership and commitment of these rows belong to the nonlinear state
 * manager rather than to the constitutive model.
 *
 * @see material::Material
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../core/core.h"

#include <memory>

namespace fem {

struct AxialStrainLinearized;
struct AxialStrainGreenLagrange;
struct VolumeStrainLinearized;
struct VolumeStrainGreenLagrange;
struct BeamGeneralizedStrain;
struct ShellMaterialStrainLinearized;
struct ShellMaterialStrainGreenLagrange;

struct AxialStressCauchy;
struct AxialStressPK2;
struct VolumeStressCauchy;
struct VolumeStressPK2;
struct BeamStressResultants;
struct ShellMaterialStressCauchy;
struct ShellMaterialStressPK2;

namespace material {

/**
 * @brief Polymorphic interface for axial, solid, beam and shell elasticity.
 *
 * The interface exposes separate overloads for linearized and finite-strain
 * evaluations so their stress measures remain explicit. Linearized axial,
 * volume and shell calls return Cauchy stress. Their finite-strain counterparts
 * accept Green-Lagrange strain and return second Piola-Kirchhoff stress. Beam
 * evaluation operates directly on generalized strains and resultants.
 *
 * Capability queries allow sections and elements to reject unsupported
 * kinematics before evaluation. The base implementations report no supported
 * formulation and raise a model error if an unsupported overload is called.
 *
 * Material-point state is passed as non-owning input/output storage.
 * `state_size()` defines the required leading components, and
 * `initialize_state()` establishes their reference history. A size of zero
 * denotes a stateless model; callers still pass valid dummy rows to keep
 * addressing uniform.
 */
struct Elasticity {
    // Shared ownership used by material definitions
    using Ptr = std::shared_ptr<Elasticity>;

    virtual ~Elasticity() = default;

    // Capability queries used by elements and sections before dispatch. Each
    // flag refers to one exact strain/stress-measure pair; support for a
    // linearized formulation does not imply support for its finite-strain form.
    virtual bool supports_axial_linearized() const;
    virtual bool supports_axial_green_lagrange() const;

    virtual bool supports_volume_linearized() const;
    virtual bool supports_volume_green_lagrange() const;

    virtual bool supports_beam_resultants() const;

    virtual bool supports_shell_integration_linearized() const;
    virtual bool supports_shell_integration_green_lagrange() const;

    // Material-point history contract. state_size() is the number of leading
    // scalar values consumed in every globally enumerated state row.
    // initialize_state() writes their constitutive reference values in place.
    virtual Index state_size() const;
    virtual void  initialize_state(Precision* state) const;

    // Every constitutive evaluation receives valid, separate state rows.
    // Implementations read old_state without modifying it and write all updated
    // history components to new_state.

    // Infinitesimal axial response in the material direction. The caller owns
    // the immutable old_state and mutable new_state rows of the current material
    // point. stress is Cauchy stress and tangent is d(sigma)/d(epsilon).
    virtual void evaluate(const AxialStrainLinearized& strain,
                          const Precision*             old_state,
                          Precision*                   new_state,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const;

    // Finite-strain axial response in the reference material direction. stress
    // is second Piola-Kirchhoff stress work-conjugate to Green-Lagrange strain;
    // tangent is the consistent derivative dS/dE.
    virtual void evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                old_state,
                          Precision*                      new_state,
                          AxialStressPK2&                 stress,
                          Precision&                      tangent) const;

    // Infinitesimal three-dimensional response in the material basis. Voigt
    // ordering follows VolumeStrain/VolumeStress, and tangent maps engineering
    // strain components to Cauchy-stress components.
    virtual void evaluate(const VolumeStrainLinearized& strain,
                          const Precision*              old_state,
                          Precision*                    new_state,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const;

    // Total-Lagrangian three-dimensional response in the reference material
    // basis. Green-Lagrange strain, PK2 stress and dS/dE remain work-conjugate;
    // the owning section handles transformations to and from global coordinates.
    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 old_state,
                          Precision*                       new_state,
                          VolumeStressPK2&                 stress,
                          Mat6&                            tangent) const;

    // Optional-tangent variant used by nonlinear residual-only assembly. A null
    // tangent requests the identical physical constitutive update and state
    // transition but permits history-dependent materials to skip an expensive
    // algorithmic tangent. The default implementation preserves compatibility by
    // evaluating the normal tangent into temporary storage.
    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 old_state,
                          Precision*                       new_state,
                          VolumeStressPK2&                 stress,
                          Mat6*                            tangent) const;

    // Generalized beam response. The section-defined six-component strain and
    // resultant ordering is preserved, and tangent is their consistent local
    // derivative. State follows the same separate input/output convention.
    virtual void evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             old_state,
                          Precision*                   new_state,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const;

    // Infinitesimal shell material response at one physical thickness point.
    // The five strain components exclude thickness-normal strain; stress is
    // Cauchy stress under the material's plane-stress reduction.
    virtual void evaluate(const ShellMaterialStrainLinearized& strain,
                          const Precision*                     old_state,
                          Precision*                           new_state,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                                 tangent) const;

    // Finite-strain shell material response at one physical thickness point.
    // The five-component Green-Lagrange input returns work-conjugate PK2 stress
    // and the consistently reduced tangent under the model's S33 = 0 convention.
    virtual void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        old_state,
                          Precision*                              new_state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5&                                   tangent) const;

    // Runtime access to concrete constitutive implementations
    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }

    template<typename T>
    const T* as() const {
        return dynamic_cast<const T*>(this);
    }
};

using ElasticityPtr = Elasticity::Ptr;

} // namespace material
} // namespace fem
