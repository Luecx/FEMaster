/**
 * @file contact.h
 * @brief Declares frictionless augmented-Lagrange contact.
 *
 * Explicit slave node regions use faceted node-to-surface contact. Slave surface
 * regions use dual-mortar segment-to-segment integration over the projected
 * slave/master overlap. Contact history is transactional so nonlinear predictor,
 * line-search and increment trials can be committed or rolled back cleanly.
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/surface.h"

#include <array>
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
 * NodeRegion slaves retain the legacy node-to-surface formulation with explicit
 * master-facet ownership. SurfaceRegion slaves use a dual-mortar formulation:
 * quadrature points only integrate the overlap, while normalized gap constraints,
 * active-set state and augmented multipliers live on global slave mortar nodes.
 */
class Contact {
    /**
     * @brief Complete nonlinear state of one contact definition.
     *
     * For NodeRegion contact `partners` stores discrete master-facet ownership.
     * For SurfaceRegion mortar contact `partners` remains empty; the current
     * coupling is reconstructed continuously from overlap integration.
     *
     * `active_slaves`, `normal_multipliers`, `gaps` and
     * `characteristic_lengths` are keyed by slave node id in both formulations.
     */
    struct AssemblyState {
        // NodeRegion only: discrete master ownership
        std::unordered_map<ID, ID> partners;

        // Unilateral state represented on slave nodes / mortar constraints
        std::unordered_set<ID>     active_slaves;
        std::unordered_map<ID, Precision> normal_multipliers;
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;

        // NodeRegion partner signature. Mortar active-set changes are handled
        // directly inside the semismooth Newton law and never set
        // last_signature_changed.
        std::uint64_t previous_signature        = 0;
        Index         previous_active           = 0;
        bool          last_signature_changed    = false;
        bool          last_augmentation_changed = false;
    };

    /**
     * @brief One nested transactional contact state.
     *
     * The freeze flags are meaningful only for NodeRegion contact. SurfaceRegion
     * mortar always recomputes overlap geometry in the current configuration.
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

    mutable RuntimeState runtime_state;

    // Low-level trial creation. Freeze policy is used only by NodeRegion contact.
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

    /**
     * Opens a trial for the next nonlinear evaluation.
     *
     * NodeRegion contact may update its master partner once and freeze it for the
     * subsequent Newton evaluations. SurfaceRegion mortar does not freeze any
     * geometric partner; only a transactional copy of the mortar state is made.
     */
    void begin_update_trial() const {
        begin_trial(false, !static_cast<bool>(slave_surfaces));
    }

    /**
     * Opens a nested trial used by predictor and line-search evaluations.
     *
     * NodeRegion partner ownership is frozen. SurfaceRegion mortar again keeps
     * only transactional state and recomputes the overlap at the trial geometry.
     */
    void begin_frozen_trial() const {
        begin_trial(!static_cast<bool>(slave_surfaces));
    }

    void commit_trial() const;
    void rollback_trial() const;

    // NodeRegion discrete partner state after a converged Newton solve
    bool partner_signature_changed() const;

    // Augmented-Lagrange updates. Surface mortar has its own update because it
    // owns no representative master partner and iterates directly over mortar
    // constraints.
    bool update_augmented_lagrange() const;
    bool update_augmented_lagrange_surface() const;

    [[nodiscard]] bool uses_slave_surface() const noexcept {
        return static_cast<bool>(slave_surfaces);
    }

    // Explicit NodeRegion slave: node-to-surface contact.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const;

    // SurfaceRegion slave: dual-mortar segment-to-segment overlap integration.
    void assemble_surface(SystemDofIds&     system_nodal_dofs,
                          model::ModelData& model_data,
                          model::NodeData&  nodal_forces,
                          TripletList&      triplets) const;
};

} // namespace constraint
} // namespace fem
