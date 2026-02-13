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
 * @brief Builds six equations to eliminate rigid-body modes for an element region.
 *
 * Accepts an element-region and internally uses all nodes of the referenced
 * elements.
 * Terms are emitted only for translational DOFs that are available in the
 * global DOF map. The equations are
 *   sum u = 0  (3 eqs)
 *   sum ( (x_i - x0) x u_i ) = 0  (3 eqs)
 * where x0 is the center of gravity of the participating structural elements.
 */
class Rbm {
public:
    // Target element region
    model::ElementRegion::Ptr element_region = nullptr;

public:
    Rbm() = default;
    explicit Rbm(model::ElementRegion::Ptr region)
        : element_region(std::move(region)) {}

    /**
     * @brief Generates 6 equations that remove rigid-body motion for the region.
     * @param system_nodal_dofs Global DOF numbering (used for availability checks).
     * @param model_data Model data providing coordinates and topology.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) const;
};

} // namespace constraint
} // namespace fem
