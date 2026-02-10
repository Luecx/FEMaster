/**
 * @file rbm.h
 * @brief Declares the RBM (rigid-body mode) constraint builder.
 *
 * The RBM constraint generates 6 well-conditioned linear equations that
 * eliminate the rigid-body motion (3 translations + 3 rotations) of a
 * component by combining selected nodal DOFs. Nodes are selected from a
 * given region and sampled to be as far apart as possible to improve
 * conditioning.
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
 * Accepts a node-region and internally collects nodes. A farthest-point
 * sampling selects up to @c max_points nodes to improve
 * conditioning. The equations are
 *   sum u = 0  (3 eqs)
 *   sum ( (x_i - x0) x u_i ) = 0  (3 eqs)
 * where x0 is the centroid of the selected nodes.
 */
class Rbm {
public:
    // Target node region
    model::NodeRegion::Ptr node_region = nullptr;

    // Sampling budget (number of nodes to use in equations)
    Index max_points = 16;

public:
    Rbm() = default;
    explicit Rbm(model::NodeRegion::Ptr region, Index max_pts = 16)
        : node_region(std::move(region)), max_points(max_pts) {}

    /**
     * @brief Generates 6 equations that remove rigid-body motion for the region.
     * @param system_nodal_dofs Global DOF numbering (used for availability checks).
     * @param model_data Model data providing coordinates and topology.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const;
};

} // namespace constraint
} // namespace fem
