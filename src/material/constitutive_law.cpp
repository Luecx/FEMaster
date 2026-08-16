/**
 * @file constitutive_law.cpp
 * @brief Implements elastic/plastic constitutive composition.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "constitutive_law.h"

#include "../core/logging.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace fem::material {

void ConstitutiveLaw::set_elasticity(ElasticityPtr elasticity) {
    logging::error(elasticity != nullptr,
        "ConstitutiveLaw requires a non-null elasticity");
    elasticity_ = std::move(elasticity);

    if (plasticity_) {
        logging::error(plasticity_->compatible_with(*elasticity_),
            "Plasticity is not compatible with the assigned elasticity");
    }
}

void ConstitutiveLaw::set_plasticity(PlasticityPtr plasticity) {
    logging::error(plasticity != nullptr,
        "ConstitutiveLaw requires a non-null plasticity");
    plasticity_ = std::move(plasticity);

    if (elasticity_) {
        logging::error(plasticity_->compatible_with(*elasticity_),
            "Plasticity is not compatible with the assigned elasticity");
    }
}

void ConstitutiveLaw::validate() const {
    logging::error(elasticity_ != nullptr,
        "ConstitutiveLaw requires elasticity");
    logging::error(!plasticity_ || plasticity_->compatible_with(*elasticity_),
        "Plasticity is not compatible with the assigned elasticity");
}

Index ConstitutiveLaw::elastic_state_size() const {
    return elasticity_ ? elasticity_->state_size() : 0;
}

Index ConstitutiveLaw::plastic_state_size() const {
    return plasticity_ ? plasticity_->state_size() : 0;
}

Index ConstitutiveLaw::state_size() const {
    return elastic_state_size() + plastic_state_size();
}

void ConstitutiveLaw::initialize_state(Precision* state) const {
    validate();
    const Index elastic_size = elastic_state_size();
    const Index plastic_size = plastic_state_size();

    logging::error(state != nullptr || elastic_size + plastic_size == 0,
        "ConstitutiveLaw received null state storage");

    if (elastic_size > 0) {
        elasticity_->initialize_state(state);
    }
    if (plastic_size > 0) {
        plasticity_->initialize_state(state + elastic_size);
    }
}

Precision* ConstitutiveLaw::prepare_elastic_trial_state(
    const Precision* state_old,
    Precision* state_new
) const {
    const Index elastic_size = elastic_state_size();
    if (elastic_size == 0) {
        return nullptr;
    }

    logging::error(state_old != nullptr && state_new != nullptr,
        "Stateful elastic backbone requires old and new state storage");
    logging::error(state_old != state_new,
        "ConstitutiveLaw state_old and state_new must not alias");

    std::copy_n(state_old, elastic_size, state_new);
    return state_new;
}

const Precision* ConstitutiveLaw::plastic_state_old(
    const Precision* state_old
) const {
    if (plastic_state_size() == 0) {
        return nullptr;
    }
    logging::error(state_old != nullptr,
        "Plastic constitutive law requires committed state storage");
    return state_old + elastic_state_size();
}

Precision* ConstitutiveLaw::plastic_state_new(Precision* state_new) const {
    if (plastic_state_size() == 0) {
        return nullptr;
    }
    logging::error(state_new != nullptr,
        "Plastic constitutive law requires trial state storage");
    return state_new + elastic_state_size();
}

bool ConstitutiveLaw::supports_axial_linearized() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_axial_linearized();
}

bool ConstitutiveLaw::supports_axial_green_lagrange() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_axial_green_lagrange();
}

bool ConstitutiveLaw::supports_volume_linearized() const {
    return elasticity_ && elasticity_->supports_volume_linearized()
        && (!plasticity_ || plasticity_->supports_volume_linearized());
}

bool ConstitutiveLaw::supports_volume_green_lagrange() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_volume_green_lagrange();
}

bool ConstitutiveLaw::supports_beam_resultants() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_beam_resultants();
}

bool ConstitutiveLaw::supports_shell_integration_linearized() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_shell_integration_linearized();
}

bool ConstitutiveLaw::supports_shell_integration_green_lagrange() const {
    return elasticity_ && !plasticity_
        && elasticity_->supports_shell_integration_green_lagrange();
}

void ConstitutiveLaw::update(const AxialStrainLinearized& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             AxialStressCauchy& stress,
                             Precision& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not support axial linearized evaluation");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::update(const AxialStrainGreenLagrange& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             AxialStressPK2& stress,
                             Precision& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not support axial Green-Lagrange evaluation");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::update(const VolumeStrainLinearized& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             VolumeStressCauchy& stress,
                             Mat6& tangent) const {
    validate();

    if (plasticity_) {
        logging::error(state_old != state_new,
            "ConstitutiveLaw state_old and state_new must not alias");
        logging::error(elastic_state_size() == 0,
            "Plasticity currently requires a stateless elastic backbone");
        plasticity_->update(
            strain,
            *elasticity_,
            plastic_state_old(state_old),
            plastic_state_new(state_new),
            stress,
            tangent
        );
        return;
    }

    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::update(const VolumeStrainGreenLagrange& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             VolumeStressPK2& stress,
                             Mat6& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Small-strain plasticity cannot be used with Green-Lagrange kinematics");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::update(const BeamGeneralizedStrain& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             BeamStressResultants& resultants,
                             Mat6& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not support generalized beam evaluation");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, resultants, tangent);
}

void ConstitutiveLaw::update(const ShellMaterialStrainLinearized& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             ShellMaterialStressCauchy& stress,
                             Mat5& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not yet support shell integration");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::update(const ShellMaterialStrainGreenLagrange& strain,
                             const Precision* state_old,
                             Precision* state_new,
                             ShellMaterialStressPK2& stress,
                             Mat5& tangent) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not yet support finite-strain shell integration");
    Precision* elastic_state = prepare_elastic_trial_state(state_old, state_new);
    elasticity_->evaluate(strain, elastic_state, stress, tangent);
}

void ConstitutiveLaw::recover(const AxialStrainLinearized& strain,
                              const Precision* state,
                              AxialStressCauchy& stress) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not support axial recovery");

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Precision tangent{};
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

void ConstitutiveLaw::recover(const AxialStrainGreenLagrange& strain,
                              const Precision* state,
                              AxialStressPK2& stress) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not support finite-strain axial recovery");

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Precision tangent{};
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

void ConstitutiveLaw::recover(const VolumeStrainLinearized& strain,
                              const Precision* state,
                              VolumeStressCauchy& stress) const {
    validate();

    if (plasticity_) {
        logging::error(elastic_state_size() == 0,
            "Plasticity currently requires a stateless elastic backbone");
        plasticity_->recover(
            strain,
            *elasticity_,
            plastic_state_old(state),
            stress
        );
        return;
    }

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Mat6 tangent;
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

void ConstitutiveLaw::recover(const VolumeStrainGreenLagrange& strain,
                              const Precision* state,
                              VolumeStressPK2& stress) const {
    validate();
    logging::error(!plasticity_,
        "Small-strain plasticity cannot recover Green-Lagrange response");

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Mat6 tangent;
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

void ConstitutiveLaw::recover(const ShellMaterialStrainLinearized& strain,
                              const Precision* state,
                              ShellMaterialStressCauchy& stress) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not yet support shell recovery");

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Mat5 tangent;
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

void ConstitutiveLaw::recover(const ShellMaterialStrainGreenLagrange& strain,
                              const Precision* state,
                              ShellMaterialStressPK2& stress) const {
    validate();
    logging::error(!plasticity_,
        "Plasticity does not yet support finite-strain shell recovery");

    const Index elastic_size = elastic_state_size();
    std::vector<Precision> scratch(static_cast<std::size_t>(elastic_size));
    if (elastic_size > 0) {
        logging::error(state != nullptr,
            "Elastic recovery requires committed state storage");
        std::copy_n(state, elastic_size, scratch.data());
    }
    Mat5 tangent;
    elasticity_->evaluate(
        strain,
        elastic_size > 0 ? scratch.data() : nullptr,
        stress,
        tangent
    );
}

} // namespace fem::material
