/**
 * @file contact.h
 * @brief Declares frictionless node-to-surface penalty contact.
 *
 * The contact definition stores one master and one slave surface region. Master
 * ownership is updated through the CalculiX-style search triangulation during a
 * Newton iteration and held fixed for all nested line-search evaluations. This
 * keeps the residual used by line search consistent with the discrete master
 * face used to build the Newton tangent.
 *
 * Contact partner state is transactional: complete increment attempts can be
 * rolled back, line-search trials can freeze the current topology, and a
 * post-Newton update can detect a discrete topology change and request a
 * discontinuity iteration at the same load factor.
 *
 * @see Contact
 * @see model::SurfaceRegion
 *
 * @author Finn Eggers
 * @date 12.08.2026
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
 * @brief Frictionless CalculiX-style node-to-surface penalty contact.
 *
 * Each unique slave node receives a representative tributary area. During an
 * update evaluation, the CalculiX-style master search graph selects one master
 * face for every contact pair retained by the positive-gap cutoff. During a
 * frozen evaluation, the stored slave-to-master connectivity is reused and only
 * the continuous closest-point geometry on that finite-element face is updated.
 *
 * The signed normal gap is
 *
 *     g = (x_s - x_m)^T n_m - clearance.
 *
 * The normal law follows the regularized linear node-to-face law used by
 * CalculiX:
 *
 *     f_n = epsilon A_s g [1/2 + atan(-g/delta)/pi].
 *
 * The slave contribution is `f_n n_m`; the opposite force is distributed to the
 * master nodes with the master shape functions at the closest point. For the
 * selected master face the tangent is analytically linearized from the
 * closest-point orthogonality equations. The discrete face ownership and
 * representative slave area are held fixed during that linearization.
 */
class Contact {
    struct TrialState {
        std::unordered_map<ID, ID> partners;
        bool allow_partner_updates = true;
        bool topology_changed      = false;
    };

    // Explicit node-to-surface contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision penalty;
    Precision clearance;
    bool      flip_normal;

    // Discrete CalculiX-style contact-element connectivity. The map contains
    // only slave nodes currently retained by the positive-gap contact cutoff.
    mutable std::unordered_map<ID, ID> partners;
    mutable bool allow_partner_updates = true;
    mutable bool topology_changed      = false;
    mutable std::vector<TrialState> trial_stack;

public:
    Contact(model::SurfaceRegion::Ptr master,
            model::SurfaceRegion::Ptr slave,
            Precision                 penalty_stiffness,
            Precision                 contact_clearance,
            bool                      flip_master_normal);

    // Assemble current node-to-surface residual and, optionally, its tangent.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList*      triplets = nullptr) const;

    // Existing model assembly passes tangent storage by reference.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const {
        assemble(system_nodal_dofs, model_data, nodal_forces, &triplets);
    }

    // Transactional discrete-contact lifecycle used by nonlinear path control.
    void begin_update_trial();
    void begin_frozen_trial();
    void commit_trial();
    void rollback_trial();

    // CalculiX freezes N2F contact-element regeneration after the early Newton
    // iterations. Continuous geometry on each retained face remains current.
    void freeze_partner_updates();

    // Reports whether the most recent update evaluation changed the discrete
    // slave-to-master contact-element topology.
    bool contact_topology_changed() const;
};

} // namespace constraint
} // namespace fem
