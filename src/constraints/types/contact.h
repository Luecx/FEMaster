/**
 * @file contact.h
 * @brief Declares frictionless dual-mortar surface-to-surface contact.
 *
 * `Contact` represents one unilateral normal-contact pair between a slave and a
 * master surface region. The physical overlap is reconstructed from the current
 * configuration during every assembly; no search radius, persistent master
 * partner, active-facet flag or quadrature-point history is stored.
 *
 * The contact law uses global slave-node augmented-Lagrange multipliers. Those
 * multipliers are transactional so nonlinear increment and line-search trials
 * can be committed or rolled back independently. The effective penalty is a
 * numerical continuation parameter: it starts from the user-provided value and
 * may increase between converged Newton solves when penetration stagnates.
 *
 * Geometric segmentation, dual-basis construction, residual/tangent assembly
 * and penalty adaptation are implemented in `contact.cpp`.
 *
 * @see Contact
 * @see model::SurfaceRegion
 * @see loadcase::tools::NonlinearStateManager
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"

#include <unordered_map>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

/**
 * @brief Frictionless dual-mortar surface-to-surface contact definition.
 *
 * Slave and master finite-element facets are projected onto one slave-centered
 * common plane, clipped there and integrated over their physical overlap. A
 * biorthogonal dual interpolation defines one normalized unilateral constraint
 * per global slave mortar node. For a nodal gap `g_i`, multiplier `lambda_i` and
 * effective penalty `epsilon`, the normal pressure coefficient is
 *
 *     p_i = max(0, lambda_i - epsilon g_i).
 *
 * The augmented-Lagrange multiplier is the only persistent physical contact
 * history. Current gaps and characteristic lengths are stored only as accepted
 * trial-evaluation data required by the subsequent post-Newton multiplier
 * update. Contact geometry itself is always recomputed from current nodal
 * positions and is therefore never committed, frozen or rolled back.
 */
class Contact {
    /**
     * @brief Transactional nodal state of one mortar contact definition.
     *
     * `multipliers` contains the persistent augmented-Lagrange normal
     * multipliers. `gaps` and `characteristic_lengths` describe the most recent
     * accepted mortar evaluation in the same nonlinear trial. They are consumed
     * by the post-Newton augmentation and may be discarded on rollback together
     * with the trial multiplier state.
     */
    struct State {
        // Persistent physical history on global slave mortar nodes
        std::unordered_map<ID, Precision> multipliers;

        // Current accepted mortar evaluation used by the next AL update
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;
    };

    // Surface-to-surface contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision         initial_penalty;
    mutable Precision penalty;
    Precision         clearance;
    bool              flip_normal;

    // Numerical penalty-adaptation state. The effective penalty persists across
    // increments and is reset only when a new nonlinear analysis is initialized.
    mutable Precision previous_augmentation_penetration = Precision(-1);
    mutable bool      previous_augmentation_changed     = false;

    // Transactional multiplier/evaluation state. Nested copies are used by
    // increment, predictor and line-search trials.
    mutable State              committed_state;
    mutable std::vector<State> trial_states;

    // Access the innermost active trial, or committed state outside a trial
    State& state() const;

public:
    // Contact construction and input parameters
    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    // Transactional nonlinear multiplier/evaluation state
    void begin_trial() const;
    void commit_trial() const;
    void rollback_trial() const;

    // Numerical augmented-Lagrange penalty adaptation
    void reset_penalty_adaptation() const;
    bool adapt_penalty() const;
    void finish_augmentation(bool changed) const;
    [[nodiscard]] Precision current_penalty() const noexcept;

    // Post-Newton augmented-Lagrange multiplier update
    bool update_augmented_lagrange() const;

    // Current mortar residual and frozen-geometry contact tangent
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
