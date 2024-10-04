#pragma once

#include "../model/sets.h"
#include "../model/surface/surface.h"

namespace fem {
namespace constraint {


class Tie {

    std::string master_set;
    std::string slave_set;

    public:
    Tie(const std::string& masterSet, const std::string& slaveSet)
        : master_set(masterSet)
        , slave_set(slaveSet) {}

    TripletList get_equations(SystemDofIds& system_nodal_dofs,
                              model::Sets<std::vector<ID>>& surface_sets,
                              model::Sets<std::vector<ID>>& node_sets,
                              std::vector<model::SurfacePtr>& surfaces,
                              NodeData& node_coords,
                              int row_offset) {
        TripletList result{};

        const std::vector<ID>& surface_ids = surface_sets.get(master_set);

        // go through each node in the slave set
        for(ID id:node_sets.get(slave_set)) {

            Vec3 node_pos;
            node_pos(0) = node_coords(id, 0);
            node_pos(1) = node_coords(id, 1);
            node_pos(2) = node_coords(id, 2);

            Precision best_dist = 1e36;
            ID best_id;
            Vec2 best_local;

            // find the closest surface
            for(ID s_id:surface_ids) {
                auto s_ptr = surfaces.at(s_id);
                if (s_ptr == nullptr) continue;

                auto local = s_ptr -> global_to_local(node_pos, node_coords);
                auto mapped = s_ptr ->local_to_global(local, node_coords);

                // compute distance
                auto dist = (node_pos - mapped).norm();

                // TODO: discard any surfaces which are too far away
                // ...

                if (dist < best_dist) {
                    best_id = s_id;
                    best_dist = dist;
                    best_local = local;
                }
            }

            auto s_ptr = surfaces.at(best_id);

            // compute which nodes to couple in the first place
            Dofs dofs_mask{};
            for(int i = 0; i < 6; i++) {
                dofs_mask(0, i) = system_nodal_dofs(id, i) >= 0;
            }
            for (int local_id = 0; local_id < s_ptr->n_nodes; local_id ++) {
                ID master_node_id = s_ptr->nodes()[local_id];
                for(ID dof_id = 0; dof_id < 6; dof_id ++) {
                    dofs_mask(0, dof_id) &= (system_nodal_dofs(master_node_id, dof_id) >= 0);
                }
            }

            auto nodal_contributions = s_ptr->shape_function(best_local);

            for(ID dof_id = 0; dof_id < 6; dof_id ++) {
                if (dofs_mask(dof_id) == false) continue;

                // set 1 for the current dof at the slave node
                result.push_back(Triplet(row_offset, system_nodal_dofs(id, dof_id), 1));

                // go through each node in the surface and couple the respective dofs
                for (int local_id = 0; local_id < s_ptr->n_nodes; local_id ++) {
                    ID master_node_id = s_ptr->nodes()[local_id];

                    Precision weight = nodal_contributions(local_id);

                    result.push_back(Triplet(row_offset, system_nodal_dofs(master_node_id, dof_id), -weight));
                }
                row_offset ++;
            }
        }

        return result;
    }

};

}
}