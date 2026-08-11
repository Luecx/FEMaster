/**
 * @file contact.h
 * @brief Declares frictionless node-to-surface penalty contact.
 *
 * The contact definition stores one master surface region and one slave surface
 * region. Slave surface nodes are treated as contact nodes; every current master
 * face is tested directly during assembly. No BVH, mortar segmentation,
 * augmented-Lagrange history, formulation switching or persistent partner state
 * is used.
 *
 * @see Contact
 * @see model::SurfaceRegion
 *
 * @author Finn Eggers
 * @date 11.08.2026
 */

#pragma once

#include "../../data/field.h"
#include "../../data/region.h"

namespace fem {
namespace model {
struct ModelData;
}

namespace constraint {

/**
 * @brief Frictionless node-to-surface penalty contact definition.
 *
 * Each unique node of the slave surface region receives a representative area
 * assembled from its incident slave facets. For every slave node, all master
 * facets are projected directly and the physically closest bounded projection is
 * selected. The signed normal gap is
 *
 *     g = (x_s - x_m)^T n_m - clearance.
 *
 * The normal law follows the regularized linear node-to-face law used by
 * CalculiX:
 *
 *     f_n = epsilon A_s g [1/2 + atan(-g/delta)/pi].
 *
 * Here `epsilon` is the penalty stiffness, `A_s` the representative slave area
 * and `delta` a small smoothing length proportional to `sqrt(A_s)`. The law is
 * evaluated up to a small positive-gap cutoff, so the transition through zero
 * gap remains smooth instead of being truncated by a hard active-set switch.
 *
 * The slave contribution is `f_n n_m`; the opposite force is distributed to the
 * master nodes with the master shape functions at the closest point. For the
 * selected master face the tangent is analytically linearized from the
 * closest-point orthogonality equations. It includes derivatives of the local
 * projection coordinates, master shape functions, surface normal, gap and smooth
 * pressure-overclosure law. The discrete master-face choice and representative
 * slave area are held fixed during one tangent evaluation, matching the
 * CalculiX node-to-face contact linearization.
 */
class Contact {
    // Explicit node-to-surface contact definition
    model::SurfaceRegion::Ptr master_surfaces;
    model::SurfaceRegion::Ptr slave_surfaces;

    Precision penalty;
    Precision clearance;
    bool      flip_normal;

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

    // Existing model assembly passes tangent storage by reference. Keep this
    // convenience overload without introducing a second contact formulation.
    void assemble(SystemDofIds&     system_nodal_dofs,
                  model::ModelData& model_data,
                  model::NodeData&  nodal_forces,
                  TripletList&      triplets) const {
        assemble(system_nodal_dofs, model_data, nodal_forces, &triplets);
    }
};

} // namespace constraint
} // namespace fem
