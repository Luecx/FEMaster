//
// Created by Finn Eggers on 04.09.23.
//
#include "../core/config.h"
#include "element/element_structural.h"
#include "model.h"
#include "shell/s8.h"

namespace fem { namespace model{
namespace {

void check_field_finite(const Field& field, const std::string& label) {
    for (Index i = 0; i < field.rows; ++i) {
        for (Index j = 0; j < field.components; ++j) {
            const bool bad = std::isnan(field(i, j)) || std::isinf(field(i, j));
            logging::error(!bad, label, " row ", i, " has invalid value at col ", j);
        }
    }
}

} // namespace

std::tuple<Field, Field>
Model::compute_stress_ip(Field& displacement, bool use_green_lagrange_nl) {
    // 1) Use IP enumeration with sentinel total in the last row
    logging::error(_data->element_ip_offsets != nullptr,
                   "element IP offset field has not been initialized");
    const auto& ip_enum = *_data->element_ip_offsets;
    logging::error(ip_enum.rows == static_cast<Index>(_data->max_elems + 1),
                   "ip_numeration must have max_elems+1 rows (with total at the end).");

    const ID total_ips = static_cast<ID>(ip_enum(static_cast<Index>(_data->max_elems), 0));
    logging::error(total_ips >= 0, "Total number of integration points must be non-negative.");

    // 2) Allocate IP-level containers (Voigt: 6 components)
    Field ip_stress{"IP_STRESS", FieldDomain::ELEMENT_IP, static_cast<Index>(total_ips), 6};
    Field ip_strain{"IP_STRAIN", FieldDomain::ELEMENT_IP, static_cast<Index>(total_ips), 6};
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

            RowMatrix rst = sel->stress_strain_ip_rst();
            logging::error(rst.rows() == sel->n_integration_points(),
                           "Element ", eid, " returned ", rst.rows(),
                           " stress IP coordinates, expected ", sel->n_integration_points());
            if (rst.rows() == 0) continue;
            sel->compute_stress_strain(&ip_strain,
                                       &ip_stress,
                                       displacement,
                                       rst,
                                       static_cast<int>(ip_offset),
                                       use_green_lagrange_nl);
        }
    }

    check_field_finite(ip_stress, "IP stress");
    check_field_finite(ip_strain, "IP strain");

    return {ip_stress, ip_strain};
}

Field Model::compute_stress_state(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_ip_offsets != nullptr,
                   "element IP offset field has not been initialized");
    const auto& ip_enum = *_data->element_ip_offsets;
    logging::error(ip_enum.rows == static_cast<Index>(_data->max_elems + 1),
                   "ip_numeration must have max_elems+1 rows (with total at the end).");

    const ID total_ips = static_cast<ID>(ip_enum(static_cast<Index>(_data->max_elems), 0));
    logging::error(total_ips >= 0, "Total number of integration points must be non-negative.");

    Field ip_stress{"IP_STRESS", FieldDomain::ELEMENT_IP, static_cast<Index>(total_ips), 8};
    ip_stress.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && eid < _data->max_elems,
                           "Element id out of range in compute_stress_state: ", eid);

            const ID ip_offset = static_cast<ID>(ip_enum(static_cast<Index>(eid), 0));
            logging::error(ip_offset >= 0 && ip_offset <= total_ips,
                           "Invalid IP offset for element ", eid, ": ", ip_offset, " / total=", total_ips);

            sel->compute_stress_state(ip_stress,
                                      displacement,
                                      static_cast<int>(ip_offset),
                                      use_green_lagrange_nl);
        }
    }

    check_field_finite(ip_stress, "Stress state");

    return ip_stress;
}

std::tuple<Field, Field> Model::compute_ip_stress_strain(Field& displacement) {
    return compute_stress_ip(displacement, false);
}

Field Model::compute_ip_stress_nonlinear(Field& displacement) {
    return compute_stress_state(displacement, true);
}

Field Model::build_internal_force_nonlinear(const Field& ip_stress) {
    logging::error(ip_stress.domain == FieldDomain::ELEMENT_IP,
                   "nonlinear internal force assembly requires ELEMENT_IP stress field");
    logging::error(_data->element_ip_offsets != nullptr,
                   "element IP offset field has not been initialized");
    const auto& ip_enum = *_data->element_ip_offsets;

    Field internal{"INTERNAL_FORCES", FieldDomain::NODE, static_cast<Index>(_data->max_nodes), 6};
    internal.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && eid < _data->max_elems,
                           "Element id out of range in build_internal_force_nonlinear: ", eid);
            const ID ip_offset = static_cast<ID>(ip_enum(static_cast<Index>(eid), 0));
            sel->compute_internal_force_nonlinear(internal, ip_stress, ip_offset);
        }
    }

    for (Index i = 0; i < internal.rows; ++i) {
        for (Index j = 0; j < internal.components; ++j) {
            const bool bad = std::isnan(internal(i, j)) || std::isinf(internal(i, j));
            logging::error(!bad, "Internal force at node ", i, " has invalid value at col ", j);
        }
    }

    return internal;
}

std::tuple<Field, Field> Model::compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_nodal_offsets != nullptr,
                   "element nodal offset field has not been initialized");
    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes =
        static_cast<Index>(nodal_offsets(static_cast<Index>(_data->max_elems), 0));

    Field element_stress{"ELEMENT_NODAL_STRESS", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_strain{"ELEMENT_NODAL_STRAIN", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_ELEMENT_WEIGHTS", FieldDomain::ELEMENT, static_cast<Index>(_data->max_elems), 1};
    element_stress.set_zero();
    element_strain.set_zero();
    element_weights.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            RowMatrix rst = sel->stress_strain_nodal_rst();
            if (rst.rows() == 0) continue;
            logging::error(rst.rows() == sel->n_nodes(),
                           "Element ", sel->elem_id, " returned ", rst.rows(),
                           " nodal stress coordinates, expected ", sel->n_nodes());
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_stress_strain(&element_strain,
                                       &element_stress,
                                       displacement,
                                       rst,
                                       static_cast<int>(offset),
                                       use_green_lagrange_nl);
            element_weights(static_cast<Index>(sel->elem_id), 0) = Precision(1);
        }
    }

    Field stress = _data->element_nodal_to_nodal(element_stress, element_weights, "STRESS");
    Field strain = _data->element_nodal_to_nodal(element_strain, element_weights, "STRAIN");
    check_field_finite(stress, "Nodal stress");
    check_field_finite(strain, "Nodal strain");

    return {stress, strain};
}

std::tuple<Field, Field> Model::compute_stress_strain(Field& displacement) {
    return compute_stress_nodal(displacement, false);
}

std::tuple<Field, Field> Model::compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_nodal_offsets != nullptr,
                   "element nodal offset field has not been initialized");
    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes =
        static_cast<Index>(nodal_offsets(static_cast<Index>(_data->max_elems), 0));
    const Index total_elements      = static_cast<Index>(_data->max_elems);

    Field element_top    {"ELEMENT_NODAL_STRESS_TOP"      , FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_bot    {"ELEMENT_NODAL_STRESS_BOT"      , FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_TOP_BOT_ELEMENT_WEIGHTS", FieldDomain::ELEMENT      , total_elements     , 1};
    element_top.set_zero();
    element_bot.set_zero();
    element_weights.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            // r,s,t values of nodes of the element
            RowMatrix base_rst = sel->stress_strain_nodal_rst();
            if (base_rst.rows() == 0) continue;

            const bool is_shell = sel->is_shell();
            RowMatrix rst_bot = base_rst;
            RowMatrix rst_top = base_rst;
            if (is_shell) {
                for (int i = 0; i < base_rst.rows(); ++i) {
                    rst_bot(i, 2) = -1;
                    rst_top(i, 2) =  1;
                }
            }

            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_stress_strain(nullptr, &element_bot, displacement, rst_bot, static_cast<int>(offset), use_green_lagrange_nl);
            sel->compute_stress_strain(nullptr, &element_top, displacement, rst_top, static_cast<int>(offset), use_green_lagrange_nl);
            element_weights(static_cast<Index>(sel->elem_id), 0) = Precision(1);
        }
    }

    Field stress_top = _data->element_nodal_to_nodal(element_top, element_weights, "STRESS_TOP");
    Field stress_bot = _data->element_nodal_to_nodal(element_bot, element_weights, "STRESS_BOT");
    check_field_finite(stress_top, "Nodal top stress");
    check_field_finite(stress_bot, "Nodal bottom stress");

    return {stress_top, stress_bot};
}

std::tuple<Field, Field> Model::compute_shell_stress_surfaces(Field& displacement) {
    return compute_stress_top_bot(displacement, false);
}

Field Model::compute_shell_resultants(Field& displacement) {
    Field resultants{"SHELL_RESULTANTS", FieldDomain::NODE, static_cast<Index>(_data->max_nodes), 8};
    resultants.set_zero();

    Field count{"SHELL_RESULTANTS_COUNT", FieldDomain::NODE, static_cast<Index>(_data->max_nodes), 1};
    count.set_zero();

    for (auto el: _data->elements) {
        if (el == nullptr) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_shell_section_forces(resultants, count, displacement);
        }
    }

    const Index max_nodes = static_cast<Index>(_data->max_nodes);

    for (Index i = 0; i < max_nodes; i++) {
        if (count(i, 0) != Precision(0)) {
            for (Index j = 0; j < resultants.components; j++) {
                resultants(i, j) /= count(i, 0);
            }
        }
    }

    for (Index i = 0; i < max_nodes; i++) {
        for (Index j = 0; j < resultants.components; j++) {
            const bool invalid = std::isnan(resultants(i, j)) || std::isinf(resultants(i, j));
            logging::error(!invalid, "Node ", i, " has nan or inf shell resultant at col ", j,
                           ". Node Usage=", count(i, 0));
        }
    }

    return resultants;
}

Field Model::compute_compliance(Field& displacement) {
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

Field Model::compute_compliance_angle_derivative(Field& displacement) {
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
Field
Model::compute_section_forces(Field& displacement) {
    logging::error(_data->element_nodal_offsets != nullptr,
                   "element nodal offset field has not been initialized");
    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes =
        static_cast<Index>(nodal_offsets(static_cast<Index>(_data->max_elems), 0));

    Field beam_forces{"BEAM_SECTION_FORCES", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    beam_forces.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_beam_section_forces(beam_forces, displacement, static_cast<int>(offset));
        }
    }

    return beam_forces;
}

Field
Model::compute_shear_flow(Field& displacement) {
    logging::error(_data->element_nodal_offsets != nullptr,
                   "element nodal offset field has not been initialized");
    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes =
        static_cast<Index>(nodal_offsets(static_cast<Index>(_data->max_elems), 0));

    Field shear_flow{"SHEAR_FLOW", FieldDomain::ELEMENT_NODAL, total_element_nodes, 1};
    shear_flow.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_shear_flow(shear_flow, displacement, static_cast<int>(offset));
        }
    }

    return shear_flow;
}
} }
