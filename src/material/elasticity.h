/**
 * @file elasticity.h
 * @brief Declares the common interface for elastic constitutive models.
 *
 * `Elasticity` separates element and section kinematics from constitutive
 * evaluation. Callers select the overload matching their strain measure and
 * dimensional reduction; implementations return the work-conjugate stress and
 * consistent tangent in the same material basis.
 *
 * Every evaluation receives the mutable state row of the current material
 * point. Stateless elastic laws leave that row unchanged, while history-
 * dependent implementations may update it in place. Ownership, reset and
 * commitment of these rows belong to the nonlinear state manager rather than
 * to the constitutive model.
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
 * Material-point state is passed as mutable, non-owning storage. `state_size()`
 * defines the required leading components, and `initialize_state()` establishes
 * their reference history. A size of zero denotes a stateless model; callers
 * may still pass a valid dummy row to keep addressing uniform.
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

    // Infinitesimal axial response in the material direction. The caller owns
    // state and passes the row of the current material point. stress is Cauchy
    // stress and tangent is d(sigma)/d(epsilon).
    virtual void evaluate(const AxialStrainLinearized& strain,
                          Precision*                   state,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const;

    // Finite-strain axial response in the reference material direction. stress
    // is second Piola-Kirchhoff stress work-conjugate to Green-Lagrange strain;
    // tangent is the consistent derivative dS/dE.
    virtual void evaluate(const AxialStrainGreenLagrange& strain,
                          Precision*                      state,
                          AxialStressPK2&                 stress,
                          Precision&                      tangent) const;

    // Infinitesimal three-dimensional response in the material basis. Voigt
    // ordering follows VolumeStrain/VolumeStress, and tangent maps engineering
    // strain components to Cauchy-stress components.
    virtual void evaluate(const VolumeStrainLinearized& strain,
                          Precision*                    state,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const;

    // Total-Lagrangian three-dimensional response in the reference material
    // basis. Green-Lagrange strain, PK2 stress and dS/dE remain work-conjugate;
    // the owning section handles transformations to and from global coordinates.
    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          Precision*                       state,
                          VolumeStressPK2&                 stress,
                          Mat6&                            tangent) const;

    // Generalized beam response. The section-defined six-component strain and
    // resultant ordering is preserved, and tangent is their consistent local
    // derivative. The state row has the same non-owning in-place semantics.
    virtual void evaluate(const BeamGeneralizedStrain& strain,
                          Precision*                   state,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const;

    // Infinitesimal shell material response at one physical thickness point.
    // The five strain components exclude thickness-normal strain; stress is
    // Cauchy stress under the material's plane-stress reduction.
    virtual void evaluate(const ShellMaterialStrainLinearized& strain,
                          Precision*                           state,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                                 tangent) const;

    // Finite-strain shell material response at one physical thickness point.
    // The five-component Green-Lagrange input returns work-conjugate PK2 stress
    // and the consistently reduced tangent under the model's S33 = 0 convention.
    virtual void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          Precision*                              state,
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
