/**
 * @file plasticity.h
 * @brief Declares history-dependent plastic constitutive components.
 *
 * Plasticity is a component of a ConstitutiveLaw. It owns the material-point
 * history contract and receives a stateless elastic law for the elastic
 * predictor. Trial updates read an immutable committed state and write a
 * separate trial state; recovery is strictly state-neutral.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "../core/core.h"

#include <memory>

namespace fem {

struct VolumeStrainLinearized;
struct VolumeStressCauchy;

namespace material {

struct Elasticity;

struct Plasticity {
    using Ptr = std::shared_ptr<Plasticity>;

    virtual ~Plasticity() = default;

    virtual Index state_size() const = 0;
    virtual void initialize_state(Precision* state) const = 0;

    virtual bool supports_volume_linearized() const { return false; }
    virtual bool compatible_with(const Elasticity& elasticity) const = 0;

    virtual void update(const VolumeStrainLinearized& strain,
                        const Elasticity&              elasticity,
                        const Precision*               state_old,
                        Precision*                     state_new,
                        VolumeStressCauchy&            stress,
                        Mat6&                          tangent) const;

    virtual void recover(const VolumeStrainLinearized& strain,
                         const Elasticity&              elasticity,
                         const Precision*               state,
                         VolumeStressCauchy&            stress) const;

    template<typename T>
    T* as() { return dynamic_cast<T*>(this); }

    template<typename T>
    const T* as() const { return dynamic_cast<const T*>(this); }
};

using PlasticityPtr = Plasticity::Ptr;

} // namespace material
} // namespace fem
