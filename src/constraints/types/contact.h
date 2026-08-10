/**
 * @file contact.h
 * @brief Declares frictionless dual-mortar surface-to-surface contact.
 *
 * Contact is defined between one slave and one master surface region. The
 * current projected overlap is rebuilt during every nonlinear evaluation and
 * integrated with a dual slave interpolation. No master partner, search radius,
 * active-set flag or quadrature-point history is stored persistently.
 *
 * The unilateral normal law uses augmented-Lagrange multipliers on global slave
 * mortar nodes. Only this multiplier history and the latest evaluated nodal gaps
 * participate in nonlinear trial commit and rollback.
 *
 * @see Contact
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

/**
 * @brief Frictionless dual-mortar surface-to-surface contact.
 *
 * Slave and master finite-element facets are projected onto one slave-centered
 * common plane, clipped there and integrated over their physical overlap. A
 * biorthogonal dual basis defines one unilateral normal constraint per global
 * slave mortar node. The augmented-Lagrange multiplier is the only persistent
 * physical contact history; contact geometry is reconstructed from the current
 * nodal positions for every assembly.
 */
class Contact {
    /**
     * @brief Transactional augmented-Lagrange state of one contact definition.
     *
     * `multipliers` stores the persistent nodal normal multipliers. `gaps` and
     * `characteristic_lengths` contain the most recent mortar evaluation in the
     * same nonlinear trial and are consumed by the subsequent multiplier update.
     */
    struct State {
        std::unordered_map<ID, Precision> multipliers;
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;
    };

    // Contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::SurfaceRegion::Ptr slave_surfaces;

    mutable Precision penalty;
    Precision         clearance;
    bool              flip_normal;

    // Numerical AL penalty state. The user value is restored once at the start
    // of a nonlinear analysis. Afterwards the effective penalty persists across
    // increments and is increased only when penetration stagnates.
    mutable Precision initial_penalty                   = Precision(0);
    mutable Precision previous_augmentation_penetration = Precision(-1);
    mutable bool      previous_augmentation_changed     = false;

    // Transactional nonlinear state. Nested copies are used by increment,
    // predictor and line-search trials; geometry itself is never frozen.
    mutable State              committed_state;
    mutable std::vector<State> trial_states;

    State& state() const;

public:
    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    // Transactional nonlinear state
    void begin_trial() const;
    void commit_trial() const;
    void rollback_trial() const;

    // Penalty adaptation. The current penetration is compared with the
    // penetration after the previous multiplier augmentation. If it did not
    // decrease by at least 20%, the effective penalty is increased one decade.
    void reset_penalty_adaptation() const {
        if (!(initial_penalty > Precision(0))) {
            initial_penalty = penalty;
        }

        penalty                           = initial_penalty;
        previous_augmentation_penetration = Precision(-1);
        previous_augmentation_changed     = false;
    }

    bool adapt_penalty() const {
        State& current = state();

        Precision maximum_penetration = Precision(0);
        for (const auto& [constraint_id, gap] : current.gaps) {
            (void) constraint_id;
            if (std::isfinite(gap)) {
                maximum_penetration =
                    std::max(maximum_penetration, std::max(Precision(0), -gap));
            }
        }

        if (!previous_augmentation_changed || maximum_penetration <= Precision(0)) {
            previous_augmentation_penetration = maximum_penetration;
            return false;
        }

        const Precision previous = previous_augmentation_penetration;
        previous_augmentation_penetration = maximum_penetration;

        if (previous <= Precision(0) ||
            maximum_penetration < Precision(0.8) * previous) {
            return false;
        }

        const Precision maximum_penalty = initial_penalty * Precision(1e6);
        const Precision adapted_penalty =
            std::min(maximum_penalty, penalty * Precision(10));

        if (!(adapted_penalty > penalty) || !std::isfinite(adapted_penalty)) {
            return false;
        }

        penalty = adapted_penalty;
        return true;
    }

    void finish_augmentation(bool changed) const {
        previous_augmentation_changed = changed;
    }

    [[nodiscard]] Precision current_penalty() const noexcept {
        return penalty;
    }

    // Post-Newton augmented-Lagrange update
    bool update_augmented_lagrange() const;

    // Current surface-to-surface mortar residual and frozen-geometry tangent
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
