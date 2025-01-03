/******************************************************************************
* @file tie.h
* @brief Defines the Tie constraint class for coupling slave and master surfaces
* in FEM. The Tie class finds the closest points on master surfaces for nodes
* in a slave set and generates the equations required for coupling their DOFs.
*
* The class supports adjusting slave node positions to match master surfaces and
* couples their degrees of freedom (DOFs) based on shape functions of the master
* surface elements.
*
* @note The coupling is performed only for the DOFs that are active in both the
* slave and master nodes.
*
* @see tie.cpp
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#pragma once    // Ensures this file is only included once during compilation

#include "../data/region.h"
#include "../model/geometry/surface/surface.h"
#include "equation.h"

namespace fem {
namespace constraint {

/******************************************************************************
* @class Tie
* @brief Implements the tie constraint for coupling slave and master surfaces in FEM.
*
* The Tie class searches for the closest points on master surfaces for each node
* in the slave set and couples their DOFs. If adjustment is enabled, it modifies
* the slave node positions to match the master surface.
******************************************************************************/
class Tie {
   model::SurfaceRegion::Ptr master_set;      ///< Name of the master set (surface set)
   model::NodeRegion::Ptr slave_set;          ///< Name of the slave set (node set)
   Precision distance;      ///< Maximum allowed distance for coupling
   bool adjust;             ///< Whether to adjust slave nodes to match the master surface

   public:
   /******************************************************************************
    * @brief Constructor for the Tie class.
    *
    * Initializes the Tie object with the master and slave sets, coupling distance,
    * and whether node positions should be adjusted.
    *
    * @param masterSet Name of the master set.
    * @param slaveSet Name of the slave set.
    * @param distance Maximum allowed distance between master and slave nodes.
    * @param adjust Boolean flag indicating whether to adjust slave node positions.
    ******************************************************************************/
   Tie(const model::SurfaceRegion::Ptr masterSet,
       const model::NodeRegion::Ptr slaveSet,
       Precision distance, bool adjust);

   /******************************************************************************
    * @brief Generates the coupling equations for the tie constraint.
    *
    * This function iterates over the slave nodes, finds the closest point on the
    * master surface, and generates the coupling equations for the respective DOFs.
    * Optionally adjusts the slave node positions to match the master surface.
    *
    * @param system_nodal_dofs System-wide DOF IDs for all nodes.
    * @param model_data Model data containing node coordinates and surface geometry.
    * @return TripletList A list of triplets representing the coupling equations.
    ******************************************************************************/
   Equations get_equations(SystemDofIds& system_nodal_dofs,
                             model::ModelData& model_data);
};

}  // namespace constraint
}  // namespace fem
