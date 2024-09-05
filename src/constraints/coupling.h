//
// Created by Finn Eggers on 04.09.24.
//

#ifndef FEMASTER_COUPLING_H
#define FEMASTER_COUPLING_H

#include "../model/sets.h"

enum class CouplingType {
    KINEMATIC,
    // TODO: Add more coupling types
};

class Coupling {
    ID master_node;
    std::vector<ID> slave_nodes;
    Dof coupled_dofs;
    CouplingType type;


public:
    Coupling(ID master_node, std::vector<ID> slave_nodes, Dof coupled_dofs, CouplingType type) :
        master_node(master_node),
        slave_nodes(slave_nodes),
        coupled_dofs(coupled_dofs),
        type(type) {}

    TripletList get_equations(SystemDofIds system_nodal_dofs, int row_offset) {
        for (ID slave_node : slave_nodes) {

            for (int i = 0; i < 6; i++) {
                if (coupled_dofs(i)) {
                    system_nodal_dofs(slave_node, i) = system_nodal_dofs(master_node, i);
                }
            }
        }
    }

    DOF master_dof(SystemDof system_dof_mask) {
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

        return Dof{
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
