//
// Created by Finn Eggers on 04.09.23.
//
#include "model.h"


namespace fem::model{

std::tuple<NodeData, NodeData> Model::compute_stress_strain(NodeData& displacement){
    NodeData stress{max_nodes, 6};
    NodeData strain{max_nodes, 6};
    stress.setZero();
    strain.setZero();
    IndexVector count{max_nodes};
    count.setZero();

    for(auto el: elements){
        if(el == nullptr) continue;
        el->compute_stress(node_coords, displacement, stress, strain);
        for(int i = 0; i < el->n_nodes(); i++){
            ID id = el->nodes()[i];
            count(id) ++;
        }
    }

    for (int i = 0; i < max_nodes; i++) {
        if (count[i] != 0) {
            for (int j = 0; j < 6; j++) {
                stress(i, j) /= count[i];
                strain(i, j) /= count[i];
            }
        }
    }
    return {stress, strain};

}
ElementData Model::compute_compliance(NodeData& displacement){
    ElementData compliance{max_elements, 1};
    compliance.setZero();

    for (size_t idx = 0; idx < elements.size(); idx++) {
        auto el = elements[idx];
        if (el == nullptr) continue;
        el->compute_compliance(node_coords, displacement, compliance);
    }

    return compliance;
}

ElementData Model::compute_volumes(){

    ElementData volumes{max_elements, 1};
    volumes.setZero();

    for (size_t idx = 0; idx < elements.size(); idx++) {
        auto el = elements[idx];
        if (el == nullptr) continue;
        volumes(el->elem_id) = el->volume(node_coords);
    }

    return volumes;
}


}
