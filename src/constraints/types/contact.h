/**
 * @file contact.h
 * @brief Declares frictionless faceted node-to-surface augmented-Lagrange contact.
 *
 * Master partners and normal multipliers are stored in explicit per-contact
 * trial states. Increment, active-set, predictor and line-search evaluations can
 * therefore be committed or rolled back without hidden global history. Active
 * or augmented contact points are tracked topologically across connected master
 * facets; only new or released contact points use a global closest-point search.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/surface.h"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

/**
 * @brief Frictionless faceted node-to-surface augmented-Lagrange contact.
 *
 * The contact constraint owns all discrete runtime history required by the
 * nonlinear solver: selected master partners, active slave nodes, normal
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
     * The state combines discrete partner ownership, active slaves, augmented
     * normal multipliers and geometry recorded by the most recent assembly.
     */
    struct AssemblyState {
        // Discrete master ownership and active set
        std::unordered_map<ID, ID> partners;
        std::unordered_set<ID>     active_slaves;

        // Augmented-Lagrange state. The multiplier has the same units as
        // PENALTY * gap; slave tributary weighting is applied only to force
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
     * Geometry-independent topology and slave weights are initialized lazily and
     * reused, while `committed` and `trials` contain the nonlinear history.
     */
    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // Fixed positive slave weights initialized once from the first
        // assembled geometry.
        std::vector<ID> slave_node_ids;
        bool            slave_nodes_initialized = false;

        std::vector<std::array<ID, 4>> master_edge_neighbors;
        bool                           master_topology_initialized = false;

        std::unordered_map<ID, Precision> slave_tributary_areas;
        bool                              slave_weights_initialized = false;

        std::uint64_t call = 0;
    };

    // Contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::NodeRegion::Ptr    slave_nodes;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision distance;
    Precision penalty;
    Precision clearance;
    bool      flip_normal;

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

    // Assemble contact forces and the faceted contact tangent with fixed
    // multipliers during the inner Newton solve.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
