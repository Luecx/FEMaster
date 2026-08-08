/**
 * @file contact.h
 * @brief Declares frictionless node-to-Nagata-surface augmented-Lagrange contact.
 *
 * The contact constraint reconstructs the selected master region as one smooth
 * Nagata surface for every assembly. Closest-point locations are stored in the
 * transactional nonlinear state and followed continuously across internal
 * Nagata patches and originating finite-element surfaces.
 *
 * Contact residuals use the analytical slave gradient and locally reconstructed
 * master-position sensitivities of the projected Nagata gap. The closest-point
 * envelope property removes projection derivatives. The tangent is the
 * symmetric Gauss-Newton part of the corresponding augmented contact potential
 * and therefore omits only the gap-Hessian term.
 *
 * @see model::NagataSurface
 *
 * @author Finn Eggers
 * @date 08.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"
#include "../../model/geometry/surface/nagata_surface.h"

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
 * @brief Frictionless augmented-Lagrange contact against a smooth Nagata master.
 *
 * The selected finite-element master surfaces are reconstructed from the
 * current nodal coordinate field during every assembly. Global projection is
 * used for new or released slave nodes. Established contact points continue
 * from their stored `NagataSurface::Location` and may cross arbitrary internal
 * Nagata patch boundaries without changing physical partner ownership.
 *
 * Runtime history is transactional. Nested predictor, Newton and line-search
 * trials inherit locations, connected master components, active slaves and
 * normal multipliers and may be committed or rolled back independently.
 * Frozen trials retain both the connected master component and the active-set
 * branch so Newton and line search see one differentiable residual. A Nagata
 * patch is a local chart and is never treated as a physical contact boundary.
 *
 * The assembled residual is the gradient of the shifted quadratic normal
 * contact potential. Central differences of a locally rebuilt patch capture
 * the dependence of Nagata positions and averaged normals on the master model
 * coordinates; closest-point sensitivities follow from stationarity. The
 * symmetric tangent retains the positive-semidefinite Gauss-Newton term
 * required by the LDLT solver.
 */
class Contact {
    /**
     * @brief Complete nonlinear state of one contact configuration.
     *
     * Locations address the reconstructed Nagata topology, whose patch ordering
     * is deterministic for an unchanged master region. Components represent
     * physical connected master regions and therefore define the discrete
     * partner signature. Multipliers have pressure units before multiplication
     * by the positive slave tributary area during residual assembly.
     */
    struct AssemblyState {
        // Tracked Nagata charts and physical connected master components
        std::unordered_map<ID, model::NagataSurface::Location> locations;
        std::unordered_map<ID, model::nagata::ComponentID>      components;
        std::unordered_set<ID>                                 active_slaves;

        // Augmented-Lagrange normal multipliers without slave-area weighting
        std::unordered_map<ID, Precision> normal_multipliers;

        // Converged geometry used by the outer augmentation update
        std::unordered_map<ID, Precision> gaps;
        std::unordered_map<ID, Precision> characteristic_lengths;

        std::uint64_t previous_signature        = 0;
        Index         previous_active           = 0;
        bool          last_signature_changed    = false;
        bool          last_augmentation_changed = false;
    };

    /**
     * @brief One nested transactional contact state and its tracking policy.
     *
     * A frozen trial retains the previously selected connected master component
     * and active-set branch while still allowing continuous projection across
     * internal charts. Without an inherited location it may acquire an initial
     * chart for the state-neutral evaluation. Update trials perform one
     * component and active-set update before adopting the frozen policy.
     */
    struct TrialState {
        AssemblyState state;
        bool          freeze_master_components = false;
        bool          freeze_after_update      = false;
    };

    /**
     * @brief Persistent topology-independent data and nested nonlinear history.
     *
     * Slave membership and positive tributary weights depend only on the model
     * definition and are initialized lazily. Nagata geometry itself is rebuilt
     * from current coordinates in each assembly and is therefore not retained.
     */
    struct RuntimeState {
        AssemblyState committed;
        std::vector<TrialState> trials;

        // Fixed slave membership collected from a node or surface region
        std::vector<ID> slave_node_ids;
        bool            slave_nodes_initialized = false;

        // Positive lumped areas for a slave surface region
        std::unordered_map<ID, Precision> slave_tributary_areas;
        bool                              slave_weights_initialized = false;

        std::uint64_t call = 0;
    };

    // Contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::NodeRegion::Ptr    slave_nodes;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision penalty;
    Precision clearance;
    bool      flip_normal;

    // Persistent nonlinear runtime state
    mutable RuntimeState runtime_state;

    // Low-level trial creation used by the semantic public trial modes
    void begin_trial(bool freeze_components, bool freeze_after_update = false) const;

public:
    // Contact definitions for nodal and surface-based slave regions
    Contact(model::SurfaceRegion::Ptr master,
            model::NodeRegion::Ptr    slave,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    // Transactional contact state and connected-component tracking policy
    void begin_update_trial() const { begin_trial(false, true); }
    void begin_frozen_trial() const { begin_trial(true); }
    void commit_trial() const;
    void rollback_trial() const;

    // State changes detected after a converged inner Newton solve
    bool partner_signature_changed() const;
    bool update_augmented_lagrange() const;

    // Nagata projection, variational residual and Gauss-Newton tangent assembly
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets,
                  bool              assemble_tangent = true) const;
};

} // namespace constraint
} // namespace fem
