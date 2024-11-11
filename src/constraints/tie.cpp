/******************************************************************************
* @file tie.cpp
* @brief Implements the methods of the Tie class for applying tie constraints in FEM.
*
* This file contains the implementation of the `get_equations` method, which finds
* the closest points on master surfaces for nodes in the slave set and generates the
* coupling equations for their DOFs.
*
* @see tie.h
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#include "tie.h"

namespace fem {
namespace constraint {

/******************************************************************************
* @brief Constructor for the Tie class.
*
* Initializes the Tie object with the master and slave sets, the maximum allowed
* distance for coupling, and the flag indicating whether slave node positions
* should be adjusted.
******************************************************************************/
Tie::Tie(const model::SurfaceRegion::Ptr masterSet,
         const model::NodeRegion::Ptr slaveSet,
         Precision distance, bool adjust)
   : master_set(masterSet)
   , slave_set(slaveSet)
   , distance(distance)
   , adjust(adjust) {}

/******************************************************************************
* @brief Generates the coupling equations for the tie constraint.
*
* This function finds the closest points on the master surface for each node in
* the slave set, checks the distance to the closest point, and generates the
* corresponding coupling equations. If the `adjust` flag is set, it also adjusts
* the positions of the slave nodes to match the master surface.
*
* @param system_nodal_dofs System-wide DOF IDs for all nodes.
* @param surface_sets Set of surface IDs for the master set.
* @param node_sets Set of node IDs for the slave set.
* @param surfaces List of surfaces to search for the closest master points.
* @param node_coords Coordinates of all nodes in the system.
* @param row_offset Row offset for inserting the equations into the global system.
* @return TripletList A list of triplets representing the coupling equations.
******************************************************************************/
TripletList Tie::get_equations(SystemDofIds& system_nodal_dofs,
                              std::vector<model::SurfacePtr>& surfaces,
                              NodeData& node_coords,
                              int row_offset) {
   TripletList result{};

   // Iterate through each node in the slave set
   for(ID id : *slave_set) {
       Vec3 node_pos;
       node_pos(0) = node_coords(id, 0);
       node_pos(1) = node_coords(id, 1);
       node_pos(2) = node_coords(id, 2);

       Precision best_dist = 1e36;
       ID best_id = -1;
       Vec2 best_local;

       // Find the closest surface on the master set
       for(ID s_id : *master_set) {
           auto& s_ptr = surfaces.at(s_id);
           if (s_ptr == nullptr) continue;

           auto local = s_ptr->global_to_local(node_pos, node_coords, true);
           auto mapped = s_ptr->local_to_global(local, node_coords);

           // Compute the distance to the surface
           auto dist = (node_pos - mapped).norm();
           if (dist > distance) continue;

           if (dist < best_dist) {
               best_id = s_id;
               best_dist = dist;
               best_local = local;
           }
       }

       // Skip if no suitable surface was found
       if (best_dist > distance) continue;

       auto& s_ptr = surfaces.at(best_id);

       // Adjust slave node position if required
       if (adjust) {
           auto mapped = s_ptr->local_to_global(best_local, node_coords);
           node_coords(id, 0) = mapped(0);
           node_coords(id, 1) = mapped(1);
           node_coords(id, 2) = mapped(2);
       }

       // Compute the DOFs to be coupled
       Dofs dofs_mask{};
       for(int i = 0; i < 6; i++) {
           dofs_mask(0, i) = system_nodal_dofs(id, i) >= 0;
       }
       for (ID local_id = 0; local_id < (ID)s_ptr->n_nodes; local_id++) {
           ID master_node_id = s_ptr->nodes()[local_id];
           for(ID dof_id = 0; dof_id < 6; dof_id++) {
               dofs_mask(0, dof_id) &= (system_nodal_dofs(master_node_id, dof_id) >= 0);
           }
       }

       // Compute the nodal contributions from the master surface
       auto nodal_contributions = s_ptr->shape_function(best_local);

       // throw warning if any absolute contribution is > 10
       logging::warning(nodal_contributions.array().abs().maxCoeff() < 10 , "Nodal contributions for slave node ", id, " exceed 10.");
       logging::error  (nodal_contributions.array().abs().maxCoeff() < 100, "Nodal contributions for slave node ", id, " exceed 100.");

       // Add coupling equations
       for(ID dof_id = 0; dof_id < 6; dof_id++) {
           if (dofs_mask(dof_id) == false) continue;

           // Add the equation for the slave node
           result.push_back(Triplet(row_offset, system_nodal_dofs(id, dof_id), 1));

           // Couple with the master nodes
           for (ID local_id = 0; local_id < (ID)s_ptr->n_nodes; local_id++) {
               ID master_node_id = s_ptr->nodes()[local_id];
               Precision weight = nodal_contributions(local_id);

               result.push_back(Triplet(row_offset, system_nodal_dofs(master_node_id, dof_id), -weight));
           }

           row_offset++;
       }
   }

   return result;
}

}  // namespace constraint
}  // namespace fem
