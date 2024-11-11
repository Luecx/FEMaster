/******************************************************************************
* @file coupling.h
* @brief Defines the Coupling class for handling kinematic _couplings in FEM.
*
* The Coupling class provides functionality to couple master and slave nodes
* based on specified degrees of freedom (DOFs). It supports coupling translations
* and rotations, and generates the necessary coupling equations for FEM systems.
*
* @note Additional coupling types can be added as needed.
*
* @see coupling.cpp
* @author Finn Eggers
* @date 04.09.2024
******************************************************************************/

#pragma once  // Ensures this file is only included once during compilation

#include "../core/types.h"
#include "../data/region.h"

namespace fem {
namespace constraint {

/******************************************************************************
* @enum CouplingType
* @brief Defines the available types of _couplings in FEM.
*
* Currently, only kinematic coupling is supported, but additional types can be
* added in the future.
******************************************************************************/
enum class CouplingType {
   KINEMATIC,  ///< Kinematic coupling between master and slave nodes
   // TODO: Add more coupling types
};

/******************************************************************************
* @class Coupling
* @brief Implements a coupling constraint between master and slave nodes.
*
* The Coupling class generates the necessary equations to couple the degrees
* of freedom (DOFs) of slave nodes with a master node, based on a specified
* coupling type and the system's DOF configurations.
******************************************************************************/
class Coupling {

   public:
   ID master_node;              ///< Master node ID
   model::NodeRegion::Ptr slave_nodes; ///< List of slave node IDs
   Dofs coupled_dofs;           ///< DOFs that are to be coupled (6 DOFs per node)
   CouplingType type;           ///< Type of coupling (e.g., kinematic)

   public:
   /******************************************************************************
    * @brief Constructor for the Coupling class.
    *
    * Initializes the coupling constraint with a master node, a set of slave nodes,
    * the degrees of freedom to couple, and the type of coupling.
    *
    * @param master_node The ID of the master node.
    * @param slave_nodes A vector of IDs representing the slave nodes.
    * @param coupled_dofs A Dofs object specifying the DOFs to couple.
    * @param type The type of coupling to apply (e.g., kinematic).
    ******************************************************************************/
   Coupling(ID master_node, model::NodeRegion::Ptr slave_nodes, Dofs coupled_dofs, CouplingType type);

   /******************************************************************************
    * @brief Generates the coupling equations for the coupled DOFs.
    *
    * This function generates the necessary equations to couple the slave node
    * DOFs with the master node, based on translations and rotations. The
    * equations are returned as a list of triplets, representing matrix entries
    * in a sparse system.
    *
    * @param system_nodal_dofs The global DOF IDs for all nodes in the system.
    * @param node_coords The coordinates of all nodes in the system.
    * @param row_offset The row offset for inserting the equations into the global system.
    * @return TripletList A list of triplets representing the coupling equations.
    ******************************************************************************/
   TripletList get_equations(SystemDofIds& system_nodal_dofs, NodeData& node_coords, int row_offset);

   /******************************************************************************
    * @brief Computes the necessary DOFs for the master node based on the coupling.
    *
    * This function determines which degrees of freedom (DOFs) are required for the
    * master node, based on the coupled DOFs and the active DOFs in the slave nodes.
    *
    * @param system_dof_mask A mask of active DOFs for all nodes in the system.
    * @return Dofs A Dofs object indicating which DOFs are required for the master node.
    ******************************************************************************/
   Dofs master_dofs(SystemDofs system_dof_mask);
};

}  // namespace constraint
}  // namespace fem
