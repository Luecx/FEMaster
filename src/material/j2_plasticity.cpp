/**
 * @file j2_plasticity.cpp
 * @brief Implements small-strain associative J2 return mapping.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "j2_plasticity.h"

#include "../core/logging.h"
#include "isotropic_elasticity.h"
#include "strain/volume_strain_linearized.h"
#include "stress/volume_stress_cauchy.h"

#include <algorithm>
#include <cmath>

namespace fem::material {

void J2Plasticity::initialize_state(Precision* state) const {
    logging::error(state != nullptr,
        "J2Plasticity requires state storage");
    for (Index i = 0; i < state_components; ++i) {
        state[i] = Precision(0);
    }
}

bool J2Plasticity::compatible_with(const Elasticity& elasticity) const {
    return dynamic_cast<const IsotropicElasticity*>(&elasticity) != nullptr
        && elasticity.state_size() == 0;
}

void J2Plasticity::add_hardening_point(Precision yield_stress,
                                       Precision plastic_strain) {
    logging::error(yield_stress > Precision(0),
        "J2Plasticity requires positive yield stress");
    logging::error(plastic_strain >= Precision(0),
        "J2Plasticity requires non-negative equivalent plastic strain");

    if (hardening_.empty()) {
        logging::error(std::abs(plastic_strain) <= Precision(1e-14),
            "J2Plasticity hardening table must start at PEEQ = 0");
    } else {
        const auto& previous = hardening_.back();
        logging::error(plastic_strain > previous.plastic_strain,
            "J2Plasticity PEEQ values must be strictly increasing");
        logging::error(yield_stress >= previous.yield_stress,
            "J2Plasticity requires non-decreasing yield stress");
    }

    hardening_.push_back({yield_stress, plastic_strain});
}

void J2Plasticity::validate_hardening() const {
    logging::error(!hardening_.empty(),
        "J2Plasticity requires at least one hardening point");
}

J2Plasticity::HardeningValue J2Plasticity::hardening(Precision peeq) const {
    validate_hardening();

    if (hardening_.size() == 1) {
        return {hardening_.front().yield_stress, Precision(0)};
    }

    if (peeq >= hardening_.back().plastic_strain) {
        return {hardening_.back().yield_stress, Precision(0)};
    }

    for (Index i = 0; i + 1 < static_cast<Index>(hardening_.size()); ++i) {
        const auto& a = hardening_[static_cast<std::size_t>(i)];
        const auto& b = hardening_[static_cast<std::size_t>(i + 1)];
        if (peeq > b.plastic_strain) {
            continue;
        }

        const Precision slope = (b.yield_stress - a.yield_stress)
                              / (b.plastic_strain - a.plastic_strain);
        const Precision yield = a.yield_stress
                              + slope * (peeq - a.plastic_strain);
        return {yield, slope};
    }

    return {hardening_.back().yield_stress, Precision(0)};
}

void J2Plasticity::update(const VolumeStrainLinearized& strain,
                          const Elasticity& elasticity,
                          const Precision* state_old,
                          Precision* state_new,
                          VolumeStressCauchy& stress,
                          Mat6& tangent) const {
    validate_hardening();

    const auto* isotropic = dynamic_cast<const IsotropicElasticity*>(&elasticity);
    logging::error(isotropic != nullptr && elasticity.state_size() == 0,
        "J2Plasticity requires stateless IsotropicElasticity");
    logging::error(state_old != nullptr && state_new != nullptr,
        "J2Plasticity requires old and new state storage");
    logging::error(state_old != state_new,
        "J2Plasticity old and new state must not alias");

    Vec6 plastic_voigt_old;
    for (Index i = 0; i < plastic_strain_components; ++i) {
        plastic_voigt_old(i) = state_old[i];
    }

    const Mat3 plastic_strain_old =
        VolumeStrainLinearized(plastic_voigt_old).tensor();
    const Precision peeq_old = state_old[peeq_index];
    logging::error(peeq_old >= Precision(0),
        "J2Plasticity encountered negative committed PEEQ");

    const VolumeStrainLinearized elastic_trial_strain(
        strain.tensor() - plastic_strain_old
    );

    VolumeStressCauchy trial_stress;
    Mat6 elastic_tangent;
    elasticity.evaluate(
        elastic_trial_strain,
        nullptr,
        trial_stress,
        elastic_tangent
    );

    const Mat3 sigma_trial = trial_stress.tensor();
    const Precision pressure = sigma_trial.trace() / Precision(3);
    const Mat3 s_trial = sigma_trial - pressure * Mat3::Identity();
    const Precision q_trial = std::sqrt(
        Precision(1.5) * s_trial.squaredNorm()
    );

    const HardeningValue hardening_old = hardening(peeq_old);
    const Precision yield_function = q_trial - hardening_old.yield_stress;
    const Precision yield_tolerance = Precision(1e-12) * std::max({
        Precision(1),
        std::abs(q_trial),
        std::abs(hardening_old.yield_stress)
    });

    if (yield_function <= yield_tolerance) {
        stress = trial_stress;
        tangent = elastic_tangent;
        for (Index i = 0; i < state_components; ++i) {
            state_new[i] = state_old[i];
        }
        return;
    }

    const Precision G = isotropic->shear;
    Precision delta_peeq = yield_function
                         / (Precision(3) * G + hardening_old.slope);
    delta_peeq = std::max(delta_peeq, Precision(0));

    HardeningValue hardening_new = hardening(peeq_old + delta_peeq);
    bool local_converged = false;

    for (Index iteration = 0; iteration < 40; ++iteration) {
        hardening_new = hardening(peeq_old + delta_peeq);
        const Precision residual = q_trial
                                 - Precision(3) * G * delta_peeq
                                 - hardening_new.yield_stress;
        const Precision scale = std::max({
            Precision(1),
            std::abs(q_trial),
            std::abs(hardening_new.yield_stress)
        });

        if (std::abs(residual) <= Precision(1e-12) * scale) {
            local_converged = true;
            break;
        }

        const Precision denominator =
            Precision(3) * G + hardening_new.slope;
        logging::error(denominator > Precision(0),
            "J2Plasticity encountered non-positive return-map denominator");

        delta_peeq = std::max(
            Precision(0),
            delta_peeq + residual / denominator
        );
    }

    logging::error(local_converged,
        "J2Plasticity local return mapping did not converge");

    hardening_new = hardening(peeq_old + delta_peeq);
    logging::error(q_trial > Precision(0),
        "J2Plasticity plastic correction requires non-zero trial equivalent stress");

    const Mat3 flow_direction = Precision(1.5) * s_trial / q_trial;
    const Mat3 plastic_strain_new = plastic_strain_old
                                  + delta_peeq * flow_direction;

    const Precision beta = Precision(1)
                         - Precision(3) * G * delta_peeq / q_trial;
    const Mat3 sigma_new = pressure * Mat3::Identity() + beta * s_trial;
    stress = VolumeStressCauchy(sigma_new);

    const Vec6 plastic_voigt_new =
        VolumeStrainLinearized(plastic_strain_new).voigt();
    for (Index i = 0; i < plastic_strain_components; ++i) {
        state_new[i] = plastic_voigt_new(i);
    }
    state_new[peeq_index] = peeq_old + delta_peeq;

    // Consistent algorithmic tangent of the discrete radial-return mapping.
    // Columns are differentiated with engineering-Voigt strain perturbations so
    // the tensor/engineering shear factors remain explicit and unambiguous.
    const Precision d_delta_peeq_d_q =
        Precision(1) / (Precision(3) * G + hardening_new.slope);
    const Precision d_beta_d_q = Precision(3) * G * (
        delta_peeq / (q_trial * q_trial)
        - d_delta_peeq_d_q / q_trial
    );

    tangent.setZero();
    for (Index col = 0; col < 6; ++col) {
        Vec6 strain_direction = Vec6::Zero();
        strain_direction(col) = Precision(1);

        const Vec6 d_sigma_trial_voigt =
            elastic_tangent * strain_direction;
        const Mat3 d_sigma_trial =
            VolumeStressCauchy(d_sigma_trial_voigt).tensor();
        const Precision d_pressure = d_sigma_trial.trace() / Precision(3);
        const Mat3 d_s_trial =
            d_sigma_trial - d_pressure * Mat3::Identity();
        const Precision d_q = Precision(1.5) / q_trial
                            * (s_trial.cwiseProduct(d_s_trial)).sum();
        const Precision d_beta = d_beta_d_q * d_q;

        const Mat3 d_sigma = d_pressure * Mat3::Identity()
                           + beta * d_s_trial
                           + d_beta * s_trial;
        tangent.col(col) = VolumeStressCauchy(d_sigma).voigt();
    }

    // Associated J2 plasticity with isotropic hardening has a symmetric
    // consistent tangent. Remove only round-off asymmetry.
    tangent = Precision(0.5) * (tangent + tangent.transpose());
}

void J2Plasticity::recover(const VolumeStrainLinearized& strain,
                           const Elasticity& elasticity,
                           const Precision* state,
                           VolumeStressCauchy& stress) const {
    validate_hardening();
    logging::error(compatible_with(elasticity),
        "J2Plasticity requires stateless IsotropicElasticity");
    logging::error(state != nullptr,
        "J2Plasticity recovery requires state storage");

    Vec6 plastic_voigt;
    for (Index i = 0; i < plastic_strain_components; ++i) {
        plastic_voigt(i) = state[i];
    }

    const VolumeStrainLinearized elastic_strain(
        strain.tensor() - VolumeStrainLinearized(plastic_voigt).tensor()
    );

    Mat6 elastic_tangent;
    elasticity.evaluate(
        elastic_strain,
        nullptr,
        stress,
        elastic_tangent
    );
}

} // namespace fem::material
