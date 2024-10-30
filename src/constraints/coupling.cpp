/******************************************************************************
* @file coupling.cpp
* @brief Implements the Coupling class for handling kinematic _couplings in FEM.
*
* This file contains the implementation of the `get_equations` and `master_dofs`
* methods, which generate the coupling equations and determine the necessary DOFs
* for master nodes.
*
* @see coupling.h
* @author Finn Eggers
* @date 04.09.2024
******************************************************************************/

#include "coupling.h"

namespace fem {
namespace constraint {

/******************************************************************************
* @brief Constructor for the Coupling class.
*
* Initializes the coupling constraint with the specified master node, slave nodes,
* DOFs to couple, and the type of coupling.
******************************************************************************/
Coupling::Coupling(ID master_node, std::vector<ID> slave_nodes, Dofs coupled_dofs, CouplingType type)
   : master_node(master_node),
   slave_nodes(slave_nodes),
   coupled_dofs(coupled_dofs),
   type(type) {}

/******************************************************************************
* @brief Generates the coupling equations for the specified DOFs.
*
* This method computes the coupling equations between the master and slave nodes,
* handling both translations and rotations. The equations are returned as a list
* of triplets, which are used to populate the sparse global system matrix.
*
* For translational DOFs, the coupling is more complex and involves the rotational
* DOFs of the master node.
*
* @param system_nodal_dofs The global DOF IDs for all nodes in the system.
* @param node_coords The coordinates of all nodes in the system.
* @param row_offset The row offset for inserting the equations into the global system.
* @return TripletList A list of triplets representing the coupling equations.
******************************************************************************/
TripletList Coupling::get_equations(SystemDofIds& system_nodal_dofs, NodeData& node_coords, int row_offset) {
   TripletList triplets{};
   Index row = row_offset;

   for (ID slave_node : slave_nodes) {
       for (int i = 0; i < 6; i++) {
           if (coupled_dofs(i)) {
               if (i >= 3) {
                   // Rotational DOFs: couple with a single value
                   auto dof_master = system_nodal_dofs(master_node, i);
                   auto dof_slave = system_nodal_dofs(slave_node, i);
                   if (dof_slave < 0) continue;

                   triplets.push_back(Triplet(row, dof_slave, 1.0));
                   triplets.push_back(Triplet(row, dof_master, -1.0));
               } else {
                   // Translational DOFs: couple with master translations and rotations
                   auto dof_slave = system_nodal_dofs(slave_node, i);
                   if (dof_slave < 0) continue;

                   ID dof_master_u = system_nodal_dofs(master_node, i);
                   ID dof_master_r1 = system_nodal_dofs(master_node, (i + 1) % 3 + 3);
                   ID dof_master_r2 = system_nodal_dofs(master_node, (i + 2) % 3 + 3);

                   Precision dz = node_coords(slave_node, 2) - node_coords(master_node, 2);
                   Precision dy = node_coords(slave_node, 1) - node_coords(master_node, 1);
                   Precision dx = node_coords(slave_node, 0) - node_coords(master_node, 0);

                   Precision dr1, dr2;
                   if (i == 0) {
                       dr1 = dz;
                       dr2 = -dy;
                   } else if (i == 1) {
                       dr1 = dx;
                       dr2 = -dz;
                   } else {
                       dr1 = dy;
                       dr2 = -dx;
                   }

                   triplets.push_back(Triplet(row, dof_slave, -1.0));
                   triplets.push_back(Triplet(row, dof_master_u, 1.0));
                   triplets.push_back(Triplet(row, dof_master_r1, -dr1));
                   triplets.push_back(Triplet(row, dof_master_r2, -dr2));
               }
               row += 1;
           }
       }
   }

   return triplets;
}

/******************************************************************************
* @brief Computes the necessary DOFs for the master node based on the coupling.
*
* This method checks which degrees of freedom (DOFs) need to be active for the
* master node based on the coupling configuration and the active DOFs in the slave nodes.
*
* @param system_dof_mask A mask of active DOFs for all nodes in the system.
* @return Dofs A Dofs object indicating which DOFs are required for the master node.
******************************************************************************/
Dofs Coupling::master_dofs(SystemDofs system_dof_mask) {
   Dofs slave_dofs = {false, false, false, false, false, false};

   // Determine which DOFs are active for the slave nodes
   for (ID slave_node : slave_nodes) {
       for (int i = 0; i < 6; i++) {
           slave_dofs(i) |= system_dof_mask(slave_node, i);
       }
   }

   // Determine necessary master DOFs based on coupling and slave DOFs
   bool x_translation_needed = coupled_dofs(0) && slave_dofs(0);
   bool y_translation_needed = coupled_dofs(1) && slave_dofs(1);
   bool z_translation_needed = coupled_dofs(2) && slave_dofs(2);

   bool x_rotation_needed = (coupled_dofs(3) && slave_dofs(3)) || (slave_dofs(1) || slave_dofs(2));
   bool y_rotation_needed = (coupled_dofs(4) && slave_dofs(4)) || (slave_dofs(0) || slave_dofs(2));
   bool z_rotation_needed = (coupled_dofs(5) && slave_dofs(5)) || (slave_dofs(0) || slave_dofs(1));

   return Dofs{
       x_translation_needed,
       y_translation_needed,
       z_translation_needed,
       x_rotation_needed,
       y_rotation_needed,
       z_rotation_needed
   };
}

}  // namespace constraint
}  // namespace fem
