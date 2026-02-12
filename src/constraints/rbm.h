/**
 * @file rbm.h
 * @brief Declares the RBM (rigid-body mode) constraint builder.
 *
 * The RBM constraint generates up to 6 linear equations that eliminate the
 * rigid-body motion (3 translations + 3 rotations) of a component by combining
 * nodal DOFs from a given region.
 */

#pragma once

#include "equation.h"
#include "../data/region.h"

namespace fem {
namespace constraint {

/**
 * @class Rbm
 * @brief Builds six equations to eliminate rigid-body modes for a region.
 *
 * Accepts a node-region and internally uses all nodes of that region.
 * Terms are emitted only for translational DOFs that are available in the
 * global DOF map. The equations are
 *   sum u = 0  (3 eqs)
 *   sum ( (x_i - x0) x u_i ) = 0  (3 eqs)
 * where x0 is the centroid of the participating nodes.
 */
class Rbm {
public:
    // Target node region
    model::NodeRegion::Ptr node_region = nullptr;

public:
    Rbm() = default;
    explicit Rbm(model::NodeRegion::Ptr region)
        : node_region(std::move(region)) {}

    /**
     * @brief Generates 6 equations that remove rigid-body motion for the region.
     * @param system_nodal_dofs Global DOF numbering (used for availability checks).
     * @param model_data Model data providing coordinates and topology.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const;
};

} // namespace constraint
} // namespace fem
