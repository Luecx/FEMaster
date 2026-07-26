/**
 * @file contact.h
 * @brief Declares frictionless faceted node-to-surface penalty contact.
 *
 * Explicit slave node sets use nodal penalty stiffness. Slave surface sets are
 * reduced to unique slave nodes and use fixed positive tributary areas that are
 * initialized on the first assembly. Master-surface assignments are stored in
 * explicit per-contact trial states so line-search and increment-cutback
 * evaluations can be committed or rolled back without hidden global state.
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

class Contact {
    struct AssemblyState {
        std::unordered_map<ID, ID> partners;
        std::unordered_set<ID>     active_slaves;

        std::uint64_t previous_signature     = 0;
        Index         previous_active        = 0;
        bool          last_signature_changed = false;
    };

    struct TrialState {
        AssemblyState state;
        bool          freeze_surface_partners = false;
        bool          freeze_after_update     = false;
    };

    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // Fixed positive slave weights. They are initialized once from the
        // geometry present during the first assembly and remain constant during
        // all subsequent Newton, line-search, and cutback evaluations.
        std::vector<ID>                slave_node_ids;
        bool                           slave_nodes_initialized = false;

        std::vector<std::array<ID, 4>> master_edge_neighbors;
        bool                           master_topology_initialized = false;

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

    // Trial-state control for nonlinear solvers. Increment and active-set
    // trials may update the selected master surface once and then freeze it.
    // Natural coordinates are reprojected on that fixed surface so contact may
    // slide within the element. Line-search trials freeze immediately.
    void begin_trial(bool freeze_partners, bool freeze_after_update = false) const;
    void commit_trial() const;
    void rollback_trial() const;
    bool partner_signature_changed() const;

    // Assemble contact forces and the consistent faceted contact tangent.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem