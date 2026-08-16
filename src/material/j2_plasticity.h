/**
 * @file j2_plasticity.h
 * @brief Declares small-strain associative J2 plasticity.
 *
 * The first implementation supports three-dimensional small-strain solids with
 * IsotropicElasticity and tabulated isotropic hardening. State consists of the
 * six engineering-Voigt plastic-strain components followed by PEEQ.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "plasticity.h"

#include <vector>

namespace fem::material {

struct J2Plasticity final : Plasticity {
    struct HardeningPoint {
        Precision yield_stress = Precision(0);
        Precision plastic_strain = Precision(0);
    };

    static constexpr Index plastic_strain_components = 6;
    static constexpr Index peeq_index = 6;
    static constexpr Index state_components = 7;

    Index state_size() const override { return state_components; }
    void initialize_state(Precision* state) const override;

    bool supports_volume_linearized() const override { return true; }
    bool compatible_with(const Elasticity& elasticity) const override;

    void add_hardening_point(Precision yield_stress, Precision plastic_strain);
    const std::vector<HardeningPoint>& hardening_points() const { return hardening_; }

    void update(const VolumeStrainLinearized& strain,
                const Elasticity&              elasticity,
                const Precision*               state_old,
                Precision*                     state_new,
                VolumeStressCauchy&            stress,
                Mat6&                          tangent) const override;

    void recover(const VolumeStrainLinearized& strain,
                 const Elasticity&              elasticity,
                 const Precision*               state,
                 VolumeStressCauchy&            stress) const override;

private:
    struct HardeningValue {
        Precision yield_stress = Precision(0);
        Precision slope = Precision(0);
    };

    HardeningValue hardening(Precision peeq) const;
    void validate_hardening() const;

    std::vector<HardeningPoint> hardening_;
};

} // namespace fem::material
