/**
 * @file elasticity.h
 * @brief Declares the common interface for elastic constitutive models.
 *
 * `Elasticity` separates element and section kinematics from constitutive
 * evaluation and history integration. Read-only `evaluate()` calls inspect one
 * already established material state without changing it. `integrate()` calls
 * advance a constitutive state from an immutable source row into a separate
 * caller-owned target row.
 *
 * This source/target split is intentional: Newton iterations, line-search
 * trials, result recovery and rejected increments must never modify the last
 * accepted constitutive history implicitly.
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
 * @brief Polymorphic interface for axial, solid, beam and shell constitutive laws.
 *
 * `evaluate()` is strictly read-only with respect to material history. It is the
 * correct operation for result recovery, diagnostics and stateless response.
 * `integrate()` performs an actual constitutive update
 *
 *     (strain_{n+1}, state_n) -> (stress_{n+1}, tangent_{n+1}, state_{n+1})
 *
 * and therefore receives separate source and target state rows. Stateful models
 * must override the corresponding `integrate()` overloads. The base integration
 * implementations are valid only for stateless models and simply dispatch to
 * `evaluate()`.
 *
 * A state size of zero denotes a stateless law. Such laws accept null source and
 * target pointers. Stateful laws require valid non-aliasing source and target
 * rows whenever `integrate()` is used.
 */
struct Elasticity {
    using Ptr = std::shared_ptr<Elasticity>;

    virtual ~Elasticity() = default;

    // Capability queries used by elements and sections before dispatch.
    virtual bool supports_axial_linearized() const;
    virtual bool supports_axial_green_lagrange() const;

    virtual bool supports_volume_linearized() const;
    virtual bool supports_volume_green_lagrange() const;

    virtual bool supports_beam_resultants() const;

    virtual bool supports_shell_integration_linearized() const;
    virtual bool supports_shell_integration_green_lagrange() const;

    // Material-point history layout and initialization.
    virtual Index state_size() const;
    virtual void  initialize_state(Precision* state) const;

    // -------------------------------------------------------------------------
    // Read-only constitutive evaluation
    // -------------------------------------------------------------------------

    virtual void evaluate(const AxialStrainLinearized& strain,
                          const Precision*             state,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const;

    virtual void evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                  state,
                          AxialStressPK2&                 stress,
                          Precision&                     tangent) const;

    virtual void evaluate(const VolumeStrainLinearized& strain,
                          const Precision*              state,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const;

    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                   state,
                          VolumeStressPK2&                 stress,
                          Mat6&                           tangent) const;

    virtual void evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             state,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const;

    virtual void evaluate(const ShellMaterialStrainLinearized& strain,
                          const Precision*                     state,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                               tangent) const;

    virtual void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                          state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5&                                  tangent) const;

    // -------------------------------------------------------------------------
    // History integration from immutable source into separate target state
    // -------------------------------------------------------------------------

    virtual void integrate(const AxialStrainLinearized& strain,
                           const Precision*             state,
                           Precision*                   target_state,
                           AxialStressCauchy&           stress,
                           Precision&                   tangent) const;

    virtual void integrate(const AxialStrainGreenLagrange& strain,
                           const Precision*                  state,
                           Precision*                      target_state,
                           AxialStressPK2&                 stress,
                           Precision&                     tangent) const;

    virtual void integrate(const VolumeStrainLinearized& strain,
                           const Precision*              state,
                           Precision*                    target_state,
                           VolumeStressCauchy&           stress,
                           Mat6&                         tangent) const;

    virtual void integrate(const VolumeStrainGreenLagrange& strain,
                           const Precision*                   state,
                           Precision*                       target_state,
                           VolumeStressPK2&                 stress,
                           Mat6&                           tangent) const;

    virtual void integrate(const BeamGeneralizedStrain& strain,
                           const Precision*             state,
                           Precision*                   target_state,
                           BeamStressResultants&        resultants,
                           Mat6&                        tangent) const;

    virtual void integrate(const ShellMaterialStrainLinearized& strain,
                           const Precision*                     state,
                           Precision*                           target_state,
                           ShellMaterialStressCauchy&            stress,
                           Mat5&                               tangent) const;

    virtual void integrate(const ShellMaterialStrainGreenLagrange& strain,
                           const Precision*                          state,
                           Precision*                              target_state,
                           ShellMaterialStressPK2&                 stress,
                           Mat5&                                  tangent) const;

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
