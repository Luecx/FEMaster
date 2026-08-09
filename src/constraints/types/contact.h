/**
 * @file contact.h
 * @brief Declares frictionless augmented-Lagrange contact.
 *
 * Explicit slave node regions use node-to-surface contact. Slave surface regions
 * use segment-to-segment mortar integration over the projected overlap of slave
 * and master surface patches. Contact history is transactional so nonlinear
 * predictor, line-search and increment trials can be committed or rolled back.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/surface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

/**
 * @brief Frictionless augmented-Lagrange contact.
 *
 * NodeRegion slaves retain the faceted node-to-surface formulation. For a
 * SurfaceRegion slave the unilateral normal constraint is represented by nodal
 * mortar gaps obtained from integration over slave/master overlap segments.
 * Normal multipliers are stored per constrained slave node and interpolated over
 * the slave surface during contact-traction assembly.
 */
class Contact {
    /**
     * @brief Complete nonlinear state of one contact configuration.
     *
     * The map key is a slave node id for both explicit NodeRegion contact and
     * mortar SurfaceRegion contact. `partners` stores the selected master facet
     * for node-to-surface contact and one representative master facet for mortar;
     * mortar topology changes are detected from the unilateral nodal active set,
     * not from facet-to-facet transfer inside the overlap integration.
     */
    struct AssemblyState {
        // Discrete master ownership and unilateral active set
        std::unordered_map<ID, ID> partners;
        std::unordered_set<ID>     active_slaves;

        // Augmented normal multipliers. The multiplier has pressure-like units
        // for surface contact; physical area weighting enters only assembly.
        std::unordered_map<ID, Precision> normal_multipliers;

        // Projected normal gaps and geometric length scale used by the outer
        // augmented-Lagrange convergence/update criterion.
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;

        std::uint64_t previous_signature        = 0;
        Index         previous_active           = 0;
        bool          last_signature_changed    = false;
        bool          last_augmentation_changed = false;
    };

    /**
     * @brief One nested transactional contact state.
     */
    struct TrialState {
        AssemblyState state;
        bool          freeze_surface_partners = false;
        bool          freeze_after_update     = false;
    };

    /**
     * @brief Persistent and nested runtime data used by contact assembly.
     */
    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // NodeRegion path only: fixed slave node list and positive tributary
        // areas used by the nodal contact formulation.
        std::vector<ID> slave_node_ids;
        bool            slave_nodes_initialized = false;

        // NodeRegion path only: connected master-facet topology for partner
        // walking during sliding contact.
        std::vector<std::array<ID, 4>> master_edge_neighbors;
        bool                           master_topology_initialized = false;

        std::unordered_map<ID, Precision> slave_tributary_areas;
        bool                              slave_weights_initialized = false;

        // Solver-level numerical state for adaptive AL penalty escalation.
        Precision previous_augmentation_penetration = Precision(-1);

        std::uint64_t call = 0;
    };

    // Contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::NodeRegion::Ptr    slave_nodes;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision         distance;
    mutable Precision penalty;
    Precision         initial_penalty = penalty;
    Precision         clearance;
    bool              flip_normal;

    mutable RuntimeState runtime_state;

    // Low-level trial creation. Public callers use the semantic trial modes below.
    void begin_trial(bool freeze_partners, bool freeze_after_update = false) const;

public:
    Contact(model::SurfaceRegion::Ptr master,
            model::NodeRegion::Ptr    slave,
            Precision                 search_distance,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision                 search_distance,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    // Transactional contact state. Node-contact update trials may select a new
    // partner once and then freeze it. Mortar contact recomputes the continuous
    // overlap geometry while retaining the same transaction/rollback semantics.
    void begin_update_trial() const { begin_trial(false, true); }
    void begin_frozen_trial() const { begin_trial(true); }
    void commit_trial() const;
    void rollback_trial() const;

    // State changes detected after a converged Newton solve
    bool partner_signature_changed() const;
    bool update_augmented_lagrange() const;

    // Reset the numerical penalty before a new nonlinear solution. If the last
    // outer augmentation changed multipliers but the converged penetration did
    // not decrease by at least 20%, increase the effective penalty by one decade.
    void reset_augmented_lagrange_penalty() const {
        penalty = initial_penalty;
        runtime_state.previous_augmentation_penetration = Precision(-1);
    }

    bool adapt_augmented_lagrange_penalty() const {
        AssemblyState& state =
            runtime_state.trials.empty()
                ? runtime_state.committed
                : runtime_state.trials.back().state;

        Precision maximum_penetration = Precision(0);
        for (const auto& [node_id, gap] : state.gaps) {
            (void) node_id;
            if (std::isfinite(gap)) {
                maximum_penetration =
                    std::max(maximum_penetration, std::max(Precision(0), -gap));
            }
        }

        if (state.last_signature_changed ||
            !state.last_augmentation_changed ||
            maximum_penetration <= Precision(0)) {
            runtime_state.previous_augmentation_penetration = maximum_penetration;
            return false;
        }

        const Precision previous = runtime_state.previous_augmentation_penetration;
        runtime_state.previous_augmentation_penetration = maximum_penetration;

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

    [[nodiscard]] Precision current_penalty() const noexcept { return penalty; }

    [[nodiscard]] bool uses_slave_surface() const noexcept {
        return static_cast<bool>(slave_surfaces);
    }

    // Explicit NodeRegion slave: node-to-surface contact.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;

    // SurfaceRegion slave: segment-to-segment mortar integration over projected
    // slave/master overlap with nodal AL multipliers.
    void assemble_surface(SystemDofIds&     system_nodal_dofs,
                          model::ModelData& model_data,
                          model::NodeData&  nodal_forces,
                          TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
