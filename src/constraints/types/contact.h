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

    // Penalty continuation is reset for every new increment attempt. The user
    // value remains the starting penalty and is increased only between converged
    // Newton solves, never inside one Newton iteration.
    mutable Precision initial_penalty     = Precision(0);
    mutable Index     augmentation_count = 0;

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

    // Penalty continuation used by the nonlinear state manager
    void reset_penalty_continuation() const {
        if (!(initial_penalty > Precision(0))) {
            initial_penalty = penalty;
        }
        penalty = initial_penalty;
        augmentation_count = 0;
    }

    void advance_penalty_continuation() const {
        if (!(initial_penalty > Precision(0))) {
            initial_penalty = penalty;
        }

        ++augmentation_count;
        if (augmentation_count % 3 != 0) {
            return;
        }

        penalty = std::min(
            penalty * Precision(10),
            initial_penalty * Precision(1e4)
        );
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
