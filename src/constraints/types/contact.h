/**
 * @file contact.h
 * @brief Declares frictionless faceted augmented-Lagrange contact.
 *
 * Master partners and normal multipliers are stored in explicit per-contact
 * trial states. Increment, active-set, predictor and line-search evaluations can
 * therefore be committed or rolled back without hidden global history. Active
 * or augmented contact points are tracked topologically across connected master
 * facets; only new or released contact points use a global closest-point search.
 *
 * Explicit slave node regions retain the node-to-surface formulation. Slave
 * surface regions are evaluated at the native integration points of each slave
 * surface and their contact contributions are distributed consistently to the
 * slave nodes through the slave shape functions.
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
 * @brief Frictionless faceted augmented-Lagrange contact.
 *
 * The contact constraint owns all discrete runtime history required by the
 * nonlinear solver: selected master partners, active slave points, normal
 * multipliers and the geometry needed by the outer augmented-Lagrange update.
 * Runtime history is transactional. Nested trials allow predictor and line-search
 * evaluations to modify temporary contact state without changing their parent
 * increment state.
 *
 * A partner-update trial allows one closest-surface or topological partner update
 * and then freezes the selected master surfaces for the remaining Newton solve.
 * A frozen trial immediately inherits and freezes the current partner ownership.
 * `NonlinearStateManager` chooses between these two modes; the contact class owns
 * the actual state stack and commit/rollback semantics.
 */
class Contact {
    /**
     * @brief Complete nonlinear state of one contact configuration.
     *
     * The state combines discrete master ownership, active slave points,
     * augmented normal multipliers and geometry recorded by the most recent
     * assembly. The map key is a node id for NodeRegion slaves and a stable
     * integration-point id for SurfaceRegion slaves.
     */
    struct AssemblyState {
        // Discrete master ownership and active set
        std::unordered_map<ID, ID> partners;
        std::unordered_set<ID>     active_slaves;

        // Augmented-Lagrange state. The multiplier has the same units as
        // PENALTY * gap; slave integration weighting is applied only to force
        // assembly and is therefore not included here.
        std::unordered_map<ID, Precision> normal_multipliers;

        // Geometry recorded by the most recent assembly. These values are used
        // only by the outer augmentation update after Newton convergence.
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;

        std::uint64_t previous_signature        = 0;
        Index         previous_active           = 0;
        bool          last_signature_changed    = false;
        bool          last_augmentation_changed = false;
    };

    /**
     * @brief One nested transactional contact state.
     *
     * The trial stores a complete child state and the partner-freezing policy
     * used by subsequent contact assemblies within this transaction level.
     */
    struct TrialState {
        AssemblyState state;
        bool          freeze_surface_partners = false;
        bool          freeze_after_update     = false;
    };

    /**
     * @brief Persistent and nested runtime data used by contact assembly.
     *
     * Geometry-independent topology and node-slave weights are initialized
     * lazily and reused, while `committed` and `trials` contain nonlinear
     * history for either nodal or surface-integration-point slaves.
     */
    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // NodeRegion path only: fixed slave node list and positive tributary
        // areas used by the legacy nodal contact formulation.
        std::vector<ID> slave_node_ids;
        bool            slave_nodes_initialized = false;

        std::vector<std::array<ID, 4>> master_edge_neighbors;
        bool                           master_topology_initialized = false;

        std::unordered_map<ID, Precision> slave_tributary_areas;
        bool                              slave_weights_initialized = false;

        // Penalty adaptation is solver-level numerical state rather than contact
        // history. It is reset when a nonlinear solution starts and increased
        // only when a previous augmentation failed to reduce penetration.
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

    // Persistent nonlinear runtime state
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

    // Transactional contact state. Update trials may select partners once;
    // frozen trials preserve the current discrete partner ownership immediately.
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
    // This keeps a deliberately low user penalty usable without allowing the AL
    // loop to spend hundreds of restarts accumulating essentially the same gap.
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
        for (const auto& [point_id, gap] : state.gaps) {
            (void) point_id;
            if (std::isfinite(gap)) {
                maximum_penetration = std::max(maximum_penetration, std::max(Precision(0), -gap));
            }
        }

        if (state.last_signature_changed || !state.last_augmentation_changed || maximum_penetration <= Precision(0)) {
            runtime_state.previous_augmentation_penetration = maximum_penetration;
            return false;
        }

        const Precision previous = runtime_state.previous_augmentation_penetration;
        runtime_state.previous_augmentation_penetration = maximum_penetration;

        if (previous <= Precision(0) || maximum_penetration < Precision(0.8) * previous) {
            return false;
        }

        const Precision maximum_penalty = initial_penalty * Precision(1e6);
        const Precision adapted_penalty = std::min(maximum_penalty, penalty * Precision(10));

        if (!(adapted_penalty > penalty) || !std::isfinite(adapted_penalty)) {
            return false;
        }

        penalty = adapted_penalty;
        return true;
    }

    [[nodiscard]] Precision current_penalty() const noexcept { return penalty; }

    // Select the assembly discretisation from the resolved slave-region type.
    [[nodiscard]] bool uses_slave_surface() const noexcept {
        return static_cast<bool>(slave_surfaces);
    }

    // Explicit NodeRegion slave: existing node-to-surface formulation.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;

    // SurfaceRegion slave: evaluate contact at the surface integration points
    // and distribute residual/tangent consistently to the slave surface nodes.
    void assemble_surface(SystemDofIds&     system_nodal_dofs,
                          model::ModelData& model_data,
                          model::NodeData&  nodal_forces,
                          TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem