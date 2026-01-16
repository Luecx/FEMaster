/**
 * @file tie.h
 * @brief Declares the tie constraint relating slave nodes to master surfaces/lines.
 *
 * A tie enforces compatibility between a slave set (node set or surface set) and
 * a master geometry (surface set or line set) by projecting slave nodes onto the
 * closest master location and generating the corresponding constraint equations.
 *
 * @see src/constraints/tie.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../data/region.h"
#include "../model/geometry/line/line_interface.h"
#include "../model/geometry/surface/surface.h"
#include "equation.h"

namespace fem {
namespace constraint {

/**
 * @class Tie
 * @brief Couples slave nodes to master surfaces/lines via closest-point projection.
 *
 * Slave can be provided as:
 *  - NodeRegion (classic)
 *  - SurfaceRegion (nodes will be extracted from the surfaces)
 */
class Tie {
    // Master: either surface region or line region.
    model::SurfaceRegion::Ptr master_surfaces; ///< Master surface region (2D).
    model::LineRegion::Ptr    master_lines;    ///< Master line region (1D).

    // Slave: either node region or surface region.
    model::NodeRegion::Ptr    slave_nodes;     ///< Slave node region (direct).
    model::SurfaceRegion::Ptr slave_surfaces;  ///< Slave surface region (nodes are extracted).

    Precision distance; ///< Maximum search distance for projection.
    bool adjust;        ///< Whether to move slave nodes onto the master geometry.

public:
    // -------------------------------------------------------------------------
    // Constructors (master surfaces)
    // -------------------------------------------------------------------------

    /**
     * @brief Builds a tie constraint with a surface master and node slave.
     */
    Tie(model::SurfaceRegion::Ptr master,
        model::NodeRegion::Ptr slave,
        Precision max_distance,
        bool do_adjust);

    /**
     * @brief Builds a tie constraint with a surface master and surface slave.
     *        Slave nodes are extracted from the slave surfaces.
     */
    Tie(model::SurfaceRegion::Ptr master,
        model::SurfaceRegion::Ptr slave,
        Precision max_distance,
        bool do_adjust);

    // -------------------------------------------------------------------------
    // Constructors (master lines)
    // -------------------------------------------------------------------------

    /** Construct a tie with a line master region (1D) and node slave. */
    Tie(model::LineRegion::Ptr master,
        model::NodeRegion::Ptr slave,
        Precision max_distance,
        bool do_adjust);

    /** Construct a tie with a line master region (1D) and surface slave. */
    Tie(model::LineRegion::Ptr master,
        model::SurfaceRegion::Ptr slave,
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
