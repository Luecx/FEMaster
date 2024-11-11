//
// Created by Finn Eggers on 04.09.23.
//
#include "model.h"
#include "../core/config.h"
#include "../math/interpolate.h"

namespace fem { namespace model{

std::tuple<NodeData, NodeData> Model::compute_stress_strain(NodeData& displacement){

    NodeData stress{max_nodes, 6};
    NodeData strain{max_nodes, 6};
    stress.setZero();
    strain.setZero();
    IndexVector count{max_nodes};
    count.setZero();

    for(auto& el: elements){
        if(el == nullptr) continue;
        if(!el->is_type(StructuralType))
            continue;
        auto sel = el->as<StructuralElement>();
        sel->compute_stress_strain_nodal(node_coords, displacement, stress, strain);
        for(int i = 0; i < sel->n_nodes(); i++){
            ID id = sel->nodes()[i];
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

    // check for any nan or inf and then display the node id
    for(int i = 0; i < max_nodes; i++){
        for(int j = 0; j < 6; j++){
            bool inv_stress = std::isnan(stress(i, j)) || std::isinf(stress(i, j));
            bool inv_strain = std::isnan(strain(i, j)) || std::isinf(strain(i, j));
            logging::error(!inv_strain, "Node ", i, " has nan or inf strain. Node Usage=", count(i));
            logging::error(!inv_stress, "Node ", i, " has nan or inf stress. Node Usage=", count(i));
        }
    }

    return {stress, strain};
}

ElementData Model::compute_compliance(NodeData& displacement){
    ElementData compliance{max_elements, 1};
    compliance.setZero();

    for (size_t idx = 0; idx < elements.size(); idx++) {
        auto& el = elements[idx];
        if (el == nullptr) continue;
        if(!el->is_type(StructuralType))
            continue;
        auto sel = el->as<StructuralElement>();
        sel->compute_compliance(node_coords, displacement, compliance);
    }

    return compliance;
}

ElementData Model::compute_volumes(){

    ElementData volumes{max_elements, 1};
    volumes.setZero();

    for (size_t idx = 0; idx < elements.size(); idx++) {
        auto& el = elements[idx];
        if (el == nullptr) continue;
        if (el->is_type(StructuralType)) {
            volumes(el->elem_id) = el->as<StructuralElement>()->volume(node_coords);
        }
    }

    return volumes;
}


} }
