/**
 * @file tie.h
 * @brief Declares the tie constraint relating slave nodes to master surfaces.
 *
 * A tie enforces compatibility between a slave node set and a master surface by
 * projecting slave nodes onto the closest master location and generating the
 * corresponding constraint equations.
 *
 * @see src/constraints/tie.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../data/region.h"
#include "../model/geometry/surface/surface.h"
#include "equation.h"

namespace fem {
namespace constraint {

/**
 * @class Tie
 * @brief Couples slave nodes to master surfaces via closest-point projection.
 */
class Tie {
    // Either a surface region or a line region may be provided as master.
    model::SurfaceRegion::Ptr master_surfaces; ///< Master surface region (2D).
    model::LineRegion::Ptr    master_lines;    ///< Master line region (1D).
    model::NodeRegion::Ptr    slave_set;       ///< Slave node region.
    Precision distance;                        ///< Maximum search distance for projection.
    bool adjust;                               ///< Whether to move slave nodes onto the master geometry.

public:
    /**
     * @brief Builds a tie constraint for the given regions.
     *
     * @param master Master surface region.
     * @param slave Slave node region.
     * @param max_distance Maximum allowed projection distance.
     * @param do_adjust Whether slave nodes should be adjusted to the surface.
     */
    Tie(model::SurfaceRegion::Ptr master,
        model::NodeRegion::Ptr slave,
        Precision max_distance,
        bool do_adjust);

    /** Construct a tie with a line master region (1D). */
    Tie(model::LineRegion::Ptr master,
        model::NodeRegion::Ptr slave,
        Precision max_distance,
        bool do_adjust);

    /**
     * @brief Generates the constraint equations associated with the tie.
     *
     * @param system_nodal_dofs Global DOF numbering.
     * @param model_data Model data providing geometry information.
     * @return Equations Coupling equations enforcing the tie.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data);
};

} // namespace constraint
} // namespace fem
