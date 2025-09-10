//
// Created by Finn Eggers on 04.09.23.
//
#include "../core/config.h"
#include "../math/interpolate.h"
#include "element/element_structural.h"
#include "model.h"
#include "shell/s8.h"

namespace fem { namespace model{


std::tuple<IPData, IPData>
Model::compute_ip_stress_strain(NodeData& displacement) {
    // 1) Build IP enumeration with sentinel total in the last row
    const auto ip_enum = this->build_integration_point_numeration();
    logging::error(ip_enum.rows() == _data->max_elems + 1,
                   "ip_numeration must have max_elems+1 rows (with total at the end).");

    const ID total_ips = static_cast<ID>(ip_enum(_data->max_elems, 0));
    logging::error(total_ips >= 0, "Total number of integration points must be non-negative.");

    // 2) Allocate IP-level containers (Voigt: 6 components)
    IPData ip_stress(total_ips, 6);
    IPData ip_strain(total_ips, 6);
    ip_stress.setZero();
    ip_strain.setZero();

    // 3) Dispatch to each structural element
    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && eid < _data->max_elems,
                           "Element id out of range in compute_ip_stress_strain: ", eid);

            const ID ip_offset = static_cast<ID>(ip_enum(eid, 0));
            logging::error(ip_offset >= 0 && ip_offset <= total_ips,
                           "Invalid IP offset for element ", eid, ": ", ip_offset, " / total=", total_ips);

            sel->compute_stress_strain(ip_stress, ip_strain, displacement, ip_offset);
        }
    }

    // 4) Basic sanity checks (optional but helpful)
    for (Index i = 0; i < ip_stress.rows(); ++i) {
        for (Index j = 0; j < ip_stress.cols(); ++j) {
            const bool badS = std::isnan(ip_stress(i, j)) || std::isinf(ip_stress(i, j));
            const bool badE = std::isnan(ip_strain(i, j)) || std::isinf(ip_strain(i, j));
            logging::error(!badS, "IP ", i, " has invalid stress at col ", j);
            logging::error(!badE, "IP ", i, " has invalid strain at col ", j);
        }
    }

    return {ip_stress, ip_strain};
}

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
