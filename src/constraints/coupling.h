//
// Created by Finn Eggers on 04.09.24.
//

#ifndef FEMASTER_COUPLING_H
#define FEMASTER_COUPLING_H

#include "../model/sets.h"
#include "../core/types.h"

enum class CouplingType {
    KINEMATIC,
    // TODO: Add more coupling types
};

class Coupling {

public:
    ID master_node;
    std::vector<ID> slave_nodes;
    Dofs coupled_dofs;
    CouplingType type;


public:
    Coupling(ID master_node, std::vector<ID> slave_nodes, Dofs coupled_dofs, CouplingType type) :
        master_node(master_node),
        slave_nodes(slave_nodes),
        coupled_dofs(coupled_dofs),
        type(type) {}

    TripletList get_equations(SystemDofIds& system_nodal_dofs, NodeData& node_coords, int row_offset) {
        TripletList triplets{};
        int row = row_offset;
        for (ID slave_node : slave_nodes) {
            for (int i = 0; i < 6; i++) {
                if (coupled_dofs(i)) {

                    // if its a rotational dof, simply couple it with a single value of 1
                    if (i >= 3) {
                        auto dof_master = system_nodal_dofs(master_node, i);
                        auto dof_slave  = system_nodal_dofs(slave_node, i);
                        // if the slave does not have this dof, skip
                        if (dof_slave < 0)
                            continue;
                        triplets.push_back(Triplet(row, dof_slave ,  1.0));
                        triplets.push_back(Triplet(row, dof_master, -1.0));
                    } else {
                        auto dof_slave  = system_nodal_dofs(slave_node, i);
                        if (dof_slave < 0)
                            continue;
                        ID dof_master_u = system_nodal_dofs(master_node, i);
                        ID dof_master_r1 = system_nodal_dofs(master_node, (i+1) % 3 + 3);
                        ID dof_master_r2 = system_nodal_dofs(master_node, (i+2) % 3 + 3);

                        Precision dz = node_coords(slave_node,2) - node_coords(master_node,2);
                        Precision dy = node_coords(slave_node,1) - node_coords(master_node,1);
                        Precision dx = node_coords(slave_node,0) - node_coords(master_node,0);

                        Precision dr1;
                        Precision dr2;
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

                        triplets.push_back(Triplet(row, dof_slave   ,  1.0));
                        triplets.push_back(Triplet(row, dof_master_u, -1.0));
                        triplets.push_back(Triplet(row, dof_master_r1, dr1));
                        triplets.push_back(Triplet(row, dof_master_r2, dr2));
                    }

                    row += 1;
                }
            }
        }
        return triplets;
    }

    Dofs master_dofs(SystemDofs system_dof_mask) {
        // We will initialize the slave DOFs as all false (not active)
        Dofs slave_dofs = {false, false, false, false, false, false};

        // Check the slave nodes' DOFs to see which are active
        for (ID slave_node : slave_nodes) {
            for (int i = 0; i < 6; i++) {
                slave_dofs(i) |= system_dof_mask(slave_node, i);
            }
        }

        // Now we construct the master DOFs based on the coupling and slave DOFs
        // Translations only depend on the respective DOFs being coupled and active in the slaves
        bool x_translation_needed = coupled_dofs(0) && slave_dofs(0); // x translation
        bool y_translation_needed = coupled_dofs(1) && slave_dofs(1); // y translation
        bool z_translation_needed = coupled_dofs(2) && slave_dofs(2); // z translation

        // Rotations depend on being coupled or necessary due to translations in other axes
        bool x_rotation_needed = (coupled_dofs(3) && slave_dofs(3)) || (slave_dofs(1) || slave_dofs(2));  // x rotation needed if y or z translation is active in slaves
        bool y_rotation_needed = (coupled_dofs(4) && slave_dofs(4)) || (slave_dofs(0) || slave_dofs(2));  // y rotation needed if x or z translation is active in slaves
        bool z_rotation_needed = (coupled_dofs(5) && slave_dofs(5)) || (slave_dofs(0) || slave_dofs(1));  // z rotation needed if x or y translation is active in slaves

        return Dofs{
                x_translation_needed,  // x translation: coupled and at least one slave has x translation
                y_translation_needed,  // y translation: coupled and at least one slave has y translation
                z_translation_needed,  // z translation: coupled and at least one slave has z translation
                x_rotation_needed,     // x rotation: coupled or needed because y or z translations are active
                y_rotation_needed,     // y rotation: coupled or needed because x or z translations are active
                z_rotation_needed      // z rotation: coupled or needed because x or y translations are active
        };
    }

};


#endif //FEMASTER_COUPLING_H
