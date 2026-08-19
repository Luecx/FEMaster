//
// Created by Finn Eggers on 04.09.23.
//
#include "../core/config.h"
#include "element/element_structural.h"
#include "model.h"
#include "shell/s8.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem { namespace model {
namespace {

void check_field_finite(const Field& field, const std::string& label) {
    for (Index i = 0; i < field.rows; ++i) {
        for (Index j = 0; j < field.components; ++j) {
            const bool bad = std::isnan(field(i, j)) || std::isinf(field(i, j));
            logging::error(!bad,
                label, " row ", i, " has invalid value at col ", j);
        }
    }
}

} // namespace

Field Model::compute_stress_state(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_ip_offsets != nullptr,
        "element IP offset field has not been initialized");

    const auto& ip_enum = *_data->element_ip_offsets;
    const Index element_count = static_cast<Index>(_data->elements.size());
    const Index total_ips = _data->field_rows(FieldDomain::ELEMENT_IP);

    Field ip_stress{"IP_STRESS", FieldDomain::ELEMENT_IP, total_ips, 8};
    ip_stress.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const ID eid = sel->elem_id;
            logging::error(eid >= 0 && static_cast<Index>(eid) < element_count,
                "Element id out of range in compute_stress_state: ", eid);

            const Index ip_offset = static_cast<Index>(ip_enum(static_cast<Index>(eid), 0));
            logging::error(ip_offset <= total_ips,
                "Invalid IP offset for element ", eid, ": ", ip_offset,
                " / total=", total_ips);

            sel->compute_stress_state(
                ip_stress,
                displacement,
                static_cast<int>(ip_offset),
                use_green_lagrange_nl
            );
        }
    }

    check_field_finite(ip_stress, "Stress state");
    return ip_stress;
}

std::tuple<Field, Field> Model::compute_stress_nodal(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    Field element_stress {"ELEMENT_NODAL_STRESS", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_strain {"ELEMENT_NODAL_STRAIN", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_ELEMENT_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    element_stress.set_zero();
    element_strain.set_zero();
    element_weights.set_zero();

    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            RowMatrix rst = sel->stress_strain_nodal_rst();
            if (rst.rows() == 0) continue;

            logging::error(rst.rows() == sel->n_nodes(),
                "Element ", sel->elem_id, " returned ", rst.rows(),
                " nodal stress coordinates, expected ", sel->n_nodes());
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_stress_strain(
                &element_strain,
                &element_stress,
                displacement,
                rst,
                static_cast<int>(offset),
                use_green_lagrange_nl
            );
            element_weights(static_cast<Index>(sel->elem_id), 0) = Precision(1);
        }
    }
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);

    Field stress = _data->element_nodal_to_nodal(element_stress, element_weights, "STRESS");
    Field strain = _data->element_nodal_to_nodal(element_strain, element_weights, "STRAIN");
    check_field_finite(stress, "Nodal stress");
    check_field_finite(strain, "Nodal strain");
    return {stress, strain};
}

std::tuple<Field, Field> Model::compute_stress_top_bot(Field& displacement, bool use_green_lagrange_nl) {
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    Field element_top    {"ELEMENT_NODAL_STRESS_TOP", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_bot    {"ELEMENT_NODAL_STRESS_BOT", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    Field element_weights{"STRESS_TOP_BOT_ELEMENT_WEIGHTS", FieldDomain::ELEMENT, element_count, 1};
    element_top.set_zero();
    element_bot.set_zero();
    element_weights.set_zero();

    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;

        if (auto sel = el->as<StructuralElement>()) {
            RowMatrix base_rst = sel->stress_strain_nodal_rst();
            if (base_rst.rows() == 0) continue;

            RowMatrix rst_bot = base_rst;
            RowMatrix rst_top = base_rst;
            if (sel->is_shell()) {
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
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);

    Field stress_top = _data->element_nodal_to_nodal(element_top, element_weights, "STRESS_TOP");
    Field stress_bot = _data->element_nodal_to_nodal(element_bot, element_weights, "STRESS_BOT");
    check_field_finite(stress_top, "Nodal top stress");
    check_field_finite(stress_bot, "Nodal bottom stress");
    return {stress_top, stress_bot};
}

Field Model::compute_shell_resultants(Field& displacement) {
    const Index node_count = _data->field_rows(FieldDomain::NODE);
    Field resultants{"SHELL_RESULTANTS", FieldDomain::NODE, node_count, 8};
    Field count{"SHELL_RESULTANTS_COUNT", FieldDomain::NODE, node_count, 1};
    resultants.set_zero();
    count.set_zero();

    for (auto el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_shell_section_forces(resultants, count, displacement);
        }
    }

    for (Index i = 0; i < node_count; ++i) {
        if (count(i, 0) != Precision(0)) {
            for (Index j = 0; j < resultants.components; ++j) {
                resultants(i, j) /= count(i, 0);
            }
        }
    }

    check_field_finite(resultants, "Shell resultants");
    return resultants;
}

Field Model::compute_compliance(Field& displacement) {
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field compliance{"COMPLIANCE", FieldDomain::ELEMENT, element_count, 1};
    compliance.set_zero();

    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_compliance(displacement, compliance);
        }
    }
    return compliance;
}

Field Model::compute_compliance_angle_derivative(Field& displacement) {
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field results{"ORIENTATION_GRAD", FieldDomain::ELEMENT, element_count, 3};
    results.set_zero();

    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            sel->compute_compliance_angle_derivative(displacement, results);
        }
    }
    return results;
}

Field Model::compute_volumes() {
    const Index element_count = static_cast<Index>(_data->elements.size());
    Field volumes{"VOLUME", FieldDomain::ELEMENT, element_count, 1};
    volumes.set_zero();

    for (auto& el : _data->elements) {
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            volumes(static_cast<Index>(el->elem_id)) = sel->volume();
        }
    }
    return volumes;
}

Field Model::compute_section_forces(Field& displacement) {
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    Field beam_forces{"BEAM_SECTION_FORCES", FieldDomain::ELEMENT_NODAL, total_element_nodes, 6};
    beam_forces.set_zero();

    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_beam_section_forces(beam_forces, displacement, static_cast<int>(offset));
        }
    }
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);
    return beam_forces;
}

Field Model::compute_shear_flow(Field& displacement) {
    logging::error(_data->element_nodal_offsets != nullptr,
        "element nodal offset field has not been initialized");

    const auto& nodal_offsets = *_data->element_nodal_offsets;
    const Index total_element_nodes = _data->field_rows(FieldDomain::ELEMENT_NODAL);
    const Index element_count       = static_cast<Index>(_data->elements.size());

    Field shear_flow{"SHEAR_FLOW", FieldDomain::ELEMENT_NODAL, total_element_nodes, 1};
    shear_flow.set_zero();

    Eigen::setNbThreads(1);
#ifdef USE_MKL
    mkl_set_num_threads(1);
#endif
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, 1024) num_threads(global_config.max_threads) if(global_config.max_threads > 1)
#endif
    for (Index elem_idx = 0; elem_idx < element_count; ++elem_idx) {
        auto el = _data->elements[static_cast<std::size_t>(elem_idx)];
        if (!el) continue;
        if (auto sel = el->as<StructuralElement>()) {
            const Index offset = static_cast<Index>(nodal_offsets(static_cast<Index>(sel->elem_id), 0));
            sel->compute_shear_flow(shear_flow, displacement, static_cast<int>(offset));
        }
    }
#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
#endif
    Eigen::setNbThreads(global_config.max_threads);
    return shear_flow;
}

} }
