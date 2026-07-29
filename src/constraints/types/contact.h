/**
 * @file contact.h
 * @brief Declares frictionless node-to-Nagata-surface augmented-Lagrange contact.
 *
 * Master partners and normal multipliers are stored in explicit per-contact
 * trial states. Increment, active-set, predictor and line-search evaluations can
 * therefore be committed or rolled back without hidden global history. Active
 * or augmented contact points are tracked topologically across connected Nagata
 * patches; only new or released contact points use a global closest-point search.
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "contact_nagata.h"

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

class Contact {
    struct AssemblyState {
        // Discrete Nagata-patch ownership and active set.
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

    struct TrialState {
        AssemblyState state;
        bool          freeze_surface_partners = false;
        bool          freeze_after_update     = false;
    };

    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // Contact-only tessellation and Nagata interpolation of the master set.
        std::unique_ptr<NagataContactSurface> master_geometry;

        // Fixed positive slave weights initialized once from the first
        // assembled geometry.
        std::vector<ID> slave_node_ids;
        bool            slave_nodes_initialized = false;

        std::unordered_map<ID, Precision> slave_tributary_areas;
        bool                              slave_weights_initialized = false;

        std::uint64_t call = 0;
    };

    model::SurfaceRegion::Ptr master_surfaces;
    model::NodeRegion::Ptr    slave_nodes;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision distance;
    Precision penalty;
    Precision clearance;
    bool      flip_normal;

    mutable RuntimeState runtime_state;

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

    // Trial-state control. Line-search and predictor subtrials freeze the
    // current partner state immediately. Increment and active-set trials may
    // update partners once and freeze them afterwards.
    void begin_trial(bool freeze_partners, bool freeze_after_update = false) const;
    void commit_trial() const;
    void rollback_trial() const;

    // State changes detected after a converged Newton solve.
    bool partner_signature_changed() const;
    bool update_augmented_lagrange() const;

    // Assemble contact forces and the Nagata-patch contact tangent with fixed
    // multipliers during the inner Newton solve.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
