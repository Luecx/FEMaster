//
// Created by Finn Eggers on 04.09.23.
//
#include "../core/config.h"
#include "../math/interpolate.h"
#include "element/element_structural.h"
#include "model.h"
#include "shell/s8.h"

namespace fem { namespace model{


std::tuple<Field, Field>
Model::compute_ip_stress_strain(Field& displacement) {
    // 1) Build IP enumeration with sentinel total in the last row
    const auto ip_enum = this->build_integration_point_numeration();
    logging::error(ip_enum.rows == static_cast<Index>(_data->max_elems + 1),
                   "ip_numeration must have max_elems+1 rows (with total at the end).");

    const ID total_ips = static_cast<ID>(ip_enum(static_cast<Index>(_data->max_elems), 0));
    logging::error(total_ips >= 0, "Total number of integration points must be non-negative.");

    // 2) Allocate IP-level containers (Voigt: 6 components)
    Field ip_stress{"IP_STRESS", FieldDomain::IP, static_cast<Index>(total_ips), 6};
    Field ip_strain{"IP_STRAIN", FieldDomain::IP, static_cast<Index>(total_ips), 6};
    ip_stress.set_zero();
    ip_strain.set_zero();

    // 3) Dispatch to each structural element
    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && eid < _data->max_elems,
                           "Element id out of range in compute_ip_stress_strain: ", eid);

            const ID ip_offset = static_cast<ID>(ip_enum(static_cast<Index>(eid), 0));
            logging::error(ip_offset >= 0 && ip_offset <= total_ips,
                           "Invalid IP offset for element ", eid, ": ", ip_offset, " / total=", total_ips);

            sel->compute_stress_strain(ip_stress, ip_strain, displacement, ip_offset);
        }
    }

    // 4) Basic sanity checks (optional but helpful)
    for (Index i = 0; i < ip_stress.rows; ++i) {
        for (Index j = 0; j < ip_stress.components; ++j) {
            const bool badS = std::isnan(ip_stress(i, j)) || std::isinf(ip_stress(i, j));
            const bool badE = std::isnan(ip_strain(i, j)) || std::isinf(ip_strain(i, j));
            logging::error(!badS, "IP ", i, " has invalid stress at col ", j);
            logging::error(!badE, "IP ", i, " has invalid strain at col ", j);
        }
    }

    return {ip_stress, ip_strain};
}

std::tuple<Field, Field> Model::compute_stress_strain(Field& displacement){

    Field stress{"STRESS", FieldDomain::NODE, static_cast<Index>(_data->max_nodes), 6};
    Field strain{"STRAIN", FieldDomain::NODE, static_cast<Index>(_data->max_nodes), 6};
    stress.set_zero();
    strain.set_zero();
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

Field Model::compute_compliance(Field& displacement){
    Field compliance{"COMPLIANCE", FieldDomain::ELEMENT, static_cast<Index>(_data->max_elems), 1};
    compliance.set_zero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if(auto sel = el->as<StructuralElement>())
            sel->compute_compliance(displacement, compliance);

    }

    return compliance;
}

Field Model::compute_compliance_angle_derivative(Field& displacement){
    Field results{"ORIENTATION_GRAD", FieldDomain::ELEMENT, static_cast<Index>(_data->max_elems), 3};
    results.set_zero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if (auto sel = el->as<StructuralElement>())
            sel->compute_compliance_angle_derivative(displacement, results);
    }

    return results;
}

Field Model::compute_volumes() {
    Field volumes{"VOLUME", FieldDomain::ELEMENT, static_cast<Index>(_data->max_elems), 1};
    volumes.set_zero();

    for (size_t idx = 0; idx < _data->elements.size(); idx++) {
        auto& el = _data->elements[idx];
        if (el == nullptr) continue;
        if (auto sel = el->as<StructuralElement>())
            volumes(el->elem_id) = sel->volume();
    }

    return volumes;
}
DynamicMatrix
Model::compute_section_forces(Field& displacement) {

    using Entry = std::tuple<ID, Index, Vec6>;
    std::vector<Entry> entries;
    entries.reserve(_data->elements.size() * 2); // heuristic: mostly 2-node beams

    // 1) Collect (elem_id, local_node, Vec6) for all beam elements
    for (auto el : _data->elements) {
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {

            auto per_node_forces = sel->section_forces(displacement);

            // non-beam elements return empty vector{}; skip them
            if (per_node_forces.empty())
                continue;

            const ID eid = sel->elem_id;

            for (Index ln = 0; ln < static_cast<Index>(per_node_forces.size()); ++ln) {
                entries.emplace_back(eid, ln, per_node_forces[ln]);
            }
        }
    }

    // 2) Pack into DynamicMatrix: [elem_id, local_node, f1..f6]
    const Index n_rows = static_cast<Index>(entries.size());
    const Index n_cols = 2 + 6; // 2 index columns + 6 value columns

    DynamicMatrix mat(n_rows, n_cols);
    mat.setZero();

    for (Index i = 0; i < n_rows; ++i) {
        ID    eid;
        Index ln;
        Vec6  v;
        std::tie(eid, ln, v) = entries[static_cast<std::size_t>(i)];

        mat(i, 0) = static_cast<Precision>(eid);
        // use ln or (ln + 1) depending on whether you want 0- or 1-based in the .res file
        mat(i, 1) = static_cast<Precision>(ln);

        for (Index k = 0; k < 6; ++k) {
            mat(i, 2 + k) = v(k);
        }
    }

    return mat;
}



} }
