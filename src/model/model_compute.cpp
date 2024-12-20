//
// Created by Finn Eggers on 04.09.23.
//
#include "model.h"
#include "../core/config.h"
#include "../math/interpolate.h"

#include "element/element_structural.h"

namespace fem { namespace model{

std::tuple<NodeData, NodeData> Model::compute_stress_strain(NodeData& displacement){

    NodeData stress{_data->max_nodes, 6};
    NodeData strain{_data->max_nodes, 6};
    stress.setZero();
    strain.setZero();
    IndexVector count{_data->max_nodes};
    count.setZero();

    for(auto el: _data->elements){
        if(el == nullptr) continue;
        if(auto sel = el->as<StructuralElement>()) {
            sel->compute_stress_strain_nodal(displacement, stress, strain);
            for(int i = 0; i < sel->n_nodes(); i++){
                ID id = sel->nodes()[i];
                count(id) ++;
            }
        }
    }

    for (int i = 0; i < _data->max_nodes; i++) {
        if (count[i] != 0) {
            for (int j = 0; j < 6; j++) {
                stress(i, j) /= count[i];
                strain(i, j) /= count[i];
            }
        }
    }

    // check for any nan or inf and then display the node id
    for(int i = 0; i < _data->max_nodes; i++){
        for(int j = 0; j < 6; j++){
            bool inv_stress = std::isnan(stress(i, j)) || std::isinf(stress(i, j));
            bool inv_strain = std::isnan(strain(i, j)) || std::isinf(strain(i, j));
            logging::error(!inv_strain, "Node ", i, " has nan or inf strain. Node Usage=", count(i));
            logging::error(!inv_stress, "Node ", i, " has nan or inf stress. Node Usage=", count(i));
        }
    }

    return {stress, strain};

//    size_t n_ip = 0;
//    for(auto el: elements){
//        if(el == nullptr) continue;
//        n_ip += el->n_integration_points();
//    }
//
//    NodeData ip_stress{n_ip, 6};
//    NodeData ip_strain{n_ip, 6};
//    NodeData nd_stress{max_nodes, 6};
//    NodeData nd_strain{max_nodes, 6};
//    NodeData xyz      {n_ip, 3};
//    IndexVector n_usage{this->_data->max_nodes};
//    ip_stress.setZero();
//    ip_strain.setZero();
//    nd_stress.setZero();
//    nd_strain.setZero();
//    xyz      .setZero();
//    n_usage  .setZero();
//
//    Dim max_integration_points = 0;
//
//    n_ip = 0;
//    for(auto el: elements){
//        if(el == nullptr) continue;
//        el->compute_stress_strain(_data->get(NodeDataEntries::POSITION), displacement, ip_stress, ip_strain, xyz, n_ip);
//        for(int i = 0; i < el->n_nodes(); i++){
//            n_usage(el->nodes()[i]) ++;
//        }
//        n_ip += el->n_integration_points();
//        max_integration_points = std::max(max_integration_points, el->n_integration_points());
//    }
//
//    // matrix containing for each node the elements connecting it
//    IndexMatrix connectivity{this->_data->max_nodes, n_usage.maxCoeff()};
//    IndexVector elem_ip_ids {this->max_elements};
//    connectivity.setZero();
//    elem_ip_ids .setZero();
//
//    ID ip_id = 0;
//
//    connectivity.setConstant(-1);
//    for(auto el: elements){
//        if(el == nullptr) continue;
//        for(int i = 0; i < el->n_nodes(); i++){
//            auto n_id = el->nodes()[i];
//            connectivity(n_id, --n_usage(n_id) ) = el->elem_id;
//        }
//        elem_ip_ids(el->elem_id) = ip_id;
//        ip_id += el->n_integration_points();
//    }
//
//    // step 4: create temporary arrays used for interpolation
//    int ip_point_upperbound = connectivity.cols() * max_integration_points;
//    DynamicVector buffer{ip_point_upperbound * ip_point_upperbound};
//    buffer.setZero();
//
//#pragma omp parallel num_threads(global_config.max_threads)
//    {
//        // Thread private storage for interpolation
//        NodeData xyz_ip{ip_point_upperbound, 3};
//        NodeData val_ip{ip_point_upperbound, 12};
//        xyz_ip.setZero();
//        val_ip.setZero();
//
//#pragma omp for
//        for(int n = 0; n < max_nodes; n++){
//            // ignore unused nodes
//            if (connectivity(n,0) < 0){
//                continue;
//            }
//
//            // track the current xyz and value id
//            ID index = 0;
//            for(auto h: connectivity.row(n)){
//                if (h == -1) break;
//                auto el = elements[h];
//                auto offset = elem_ip_ids(h);
//                for(int j = 0; j < el->n_integration_points(); j++){
//                    xyz_ip.row(index) = xyz.row(offset + j);
//                    for(int m = 0; m < 6; m++){
//                        val_ip(index, m)   = ip_stress(offset + j, m);
//                        val_ip(index, m+6) = ip_strain(offset + j, m);
//                    }
//                    index++;
//                }
//            }
//            // do the actual computation
//            DynamicVector r2{12};
//            auto res = interpolator(xyz_ip.block(0,0,index, 3),  val_ip.block(0,0,index, 12), _data->get(NodeDataEntries::POSITION).row(n), &r2);
//
//            std::cout << n << " " << r2.transpose() << std::endl;
//
//            // The below memory writes should be safe as different threads write to different rows of nd_stress and nd_strain
//            nd_stress.block(n, 0, 1, 6) = res.block(0, 0, 1, 6);
//            nd_strain.block(n, 0, 1, 6) = res.block(0, 6, 1, 6);
//        }
//    } // end of omp parallel region
//
//    return {nd_stress, nd_strain};
}

ElementData Model::compute_compliance(NodeData& displacement){
    ElementData compliance{_data->max_elems, 1};
    compliance.setZero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if(auto sel = el->as<StructuralElement>())
            sel->compute_compliance(displacement, compliance);

    }

    return compliance;
}

ElementData Model::compute_compliance_angle_derivative(NodeData& displacement){
    ElementData results{_data->max_elems, 3};
    results.setZero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if (auto sel = el->as<StructuralElement>())
            sel->compute_compliance_angle_derivative(displacement, results);
    }

    return results;
}

ElementData Model::compute_volumes() {
    ElementData volumes{_data->max_elems, 1};
    volumes.setZero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if (auto sel = el->as<StructuralElement>())
            volumes(el->elem_id) = sel->volume();
    }

    return volumes;
}


} }
