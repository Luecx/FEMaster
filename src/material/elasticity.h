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

struct Elasticity {
    using Ptr = std::shared_ptr<Elasticity>;

    virtual ~Elasticity() = default;

    virtual bool supports_axial_linearized() const;
    virtual bool supports_axial_green_lagrange() const;

    virtual bool supports_volume_linearized() const;
    virtual bool supports_volume_green_lagrange() const;

    virtual bool supports_beam_resultants() const;

    virtual bool supports_shell_integration_linearized() const;
    virtual bool supports_shell_integration_green_lagrange() const;

    virtual Index state_size() const;
    virtual void  initialize_state(Precision* state) const;

    virtual void evaluate(const AxialStrainLinearized& strain,
                          const Precision*             state_old,
                          Precision*                   state_new,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const;

    virtual void evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                state_old,
                          Precision*                      state_new,
                          AxialStressPK2&                 stress,
                          Precision&                      tangent) const;

    virtual void evaluate(const VolumeStrainLinearized& strain,
                          const Precision*              state_old,
                          Precision*                    state_new,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const;

    virtual void evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 state_old,
                          Precision*                       state_new,
                          VolumeStressPK2&                 stress,
                          Mat6&                            tangent) const;

    virtual void evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             state_old,
                          Precision*                   state_new,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const;

    virtual void evaluate(const ShellMaterialStrainLinearized& strain,
                          const Precision*                     state_old,
                          Precision*                           state_new,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                                 tangent) const;

    virtual void evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        state_old,
                          Precision*                              state_new,
                          ShellMaterialStressPK2&                 stress,
                          Mat5&                                   tangent) const;

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
