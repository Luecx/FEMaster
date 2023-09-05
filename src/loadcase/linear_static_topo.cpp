//
// Created by Luecx on 04.09.2023.
//
#include "linear_static_topo.h"

fem::loadcase::LinearStaticTopo::LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model)
    : LinearStatic(id, writer, model), density(model->max_elements, 1) {
    density.setOnes();
}

void fem::loadcase::LinearStaticTopo::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    // build stiffness scalar
    auto stiffness_scalar = density.array().pow(exponent);

    auto unconstrained = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating unconstrained index matrix"
    );

    auto supp_vec = Timer::measure(
        [&]() { return this->m_model->build_support_vector(unconstrained); },
        "building support vector"
    );

    auto load_vec = Timer::measure(
        [&]() { return this->m_model->build_load_vector(unconstrained); },
        "formulating load vector"
    );

    auto stiffness = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(unconstrained, stiffness_scalar); },
        "constructing stiffness matrix"
    );

    auto impl_load_vec = Timer::measure(
        [&]() { return this->m_model->build_implicit_load_vector(stiffness, supp_vec); },
        "computing implicit load vector"
    );

    auto constrained = Timer::measure(
        [&]() { return this->m_model->build_constrained_index_matrix(supp_vec); },
        "creating constrained index matrix"
    );

    auto mapping_vec = Timer::measure(
        [&]() { return this->m_model->build_mapping_vector(unconstrained, constrained); },
        "mapping unconstrained indices to constrained ones"
    );

    load_vec += impl_load_vec;

    auto reduced_stiffness = Timer::measure(
        [&]() { return this->m_model->build_reduced_stiffness(mapping_vec, stiffness); },
        "reducing stiffness matrix"
    );

    auto reduced_load = Timer::measure(
        [&]() { return this->m_model->build_reduced_load(mapping_vec, load_vec); },
        "reducing load vector"
    );

    auto displacement = Timer::measure(
        [&]() { return solve(device, method, reduced_stiffness, reduced_load); },
        "solving system of equations"
    );

    auto disp_matrix = m_model->build_global_displacement(constrained, displacement, unconstrained, supp_vec);

    NodeData stress;
    NodeData strain;

    std::tie(stress, strain) = m_model->compute_stress_strain(disp_matrix);

    ElementData compliance_raw = m_model->compute_compliance(disp_matrix);
    ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);

    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(disp_matrix   , "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain        , "STRAIN");
    m_writer->write_eigen_matrix(stress        , "STRESS");
    m_writer->write_eigen_matrix(compliance_raw, "COMPLIANCE_RAW");
    m_writer->write_eigen_matrix(compliance_adj, "COMPLIANCE_ADJ");
    m_writer->write_eigen_matrix(dens_grad     , "DENS_GRAD");
}
