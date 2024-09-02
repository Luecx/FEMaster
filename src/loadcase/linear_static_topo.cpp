//
// Created by Luecx on 04.09.2023.
//
#include "linear_static_topo.h"

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"

#include <queue>
#include <set>

fem::loadcase::LinearStaticTopo::LinearStaticTopo(ID id, reader::Writer* writer, model::Model* model)
    : LinearStatic(id, writer, model), density(model->max_elements, 1) {
    density.setOnes();
}

void fem::loadcase::LinearStaticTopo::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC TOPO");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    auto stiffness_scalar = density.array().pow(exponent);

    // Step 1: Generate unconstrained index matrix
    auto unconstrained = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating unconstrained index matrix"
    );

    // Step 2: Build the support matrix for the constraints
    auto supp_mat = Timer::measure(
        [&]() { return this->m_model->build_constraint_matrix(supps);},
        "building support vector"
    );

    // Step 3: Build the load matrix based on applied loads
    auto load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(loads); },
        "formulating load vector"
    );

    // Step 4: Construct the global stiffness matrix
    auto stiffness = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(unconstrained, stiffness_scalar); },
        "constructing stiffness matrix"
    );

    // Step 5: Reduce load and support vectors to handle constrained degrees of freedom
    auto reduced_load_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(unconstrained, load_mat); },
        "reducing load vector"
    );

    auto reduced_supp_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(unconstrained, supp_mat); },
        "reducing support vector"
    );

    // Compute the implicit load vector to account for constraints
    auto impl_load_vec = Timer::measure(
        [&]() { return mattools::extract_scaled_row_sum(stiffness, reduced_supp_vec); },
        "computing implicit load vector"
    );

    // Add implicit load vector to the reduced load vector
    reduced_load_vec += impl_load_vec;

    // Step 6: Reduce stiffness matrix to handle constrained degrees of freedom
    auto sol_stiffness = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(stiffness, reduced_supp_vec); },
        "reducing stiffness matrix to solver"
    );

    auto sol_support = Timer::measure(
        [&]() { return mattools::reduce_vec_to_vec(reduced_supp_vec, reduced_supp_vec); },
        "reducing support vector to solver"
    );

    auto sol_load = Timer::measure(
        [&]() { return mattools::reduce_vec_to_vec(reduced_load_vec, reduced_supp_vec); },
        "reducing load vector to solver"
    );

    // Compress the stiffness matrix for efficient solving
    sol_stiffness.makeCompressed();

    // Log the overview of the system
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->max_nodes);
    logging::info(true, "total dofs        : ", unconstrained.maxCoeff() + 1);
    logging::info(true, "unconstrained dofs: ", sol_support.rows());
    logging::info(true, "constrained dofs  : ", unconstrained.maxCoeff() + 1 - sol_support.rows());
    logging::down();

    // Step 7: Solve the system of equations to get displacements
    auto sol_disp = Timer::measure(
        [&]() { return solve(device, method, sol_stiffness, sol_load); },
        "solving system of equations"
    );

    // Step 8: Expand the reduced displacement vector to full size
    auto reduced_disp = Timer::measure(
        [&]() { return mattools::expand_vec_to_vec(sol_disp, reduced_supp_vec); },
        "expanding displacement vector"
    );

    // Expand the displacement vector to match the full index matrix
    auto disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(unconstrained, reduced_disp); },
        "expanding displacement vector"
    );

    // Compute stresses and strains at the nodes
    NodeData stress;
    NodeData strain;

    std::tie(stress, strain) = Timer::measure(
        [&]() { return m_model->compute_stress_strain(disp_mat); },
        "Interpolation of stress and strain components to nodes"
    );

    ElementData compliance_raw = m_model->compute_compliance(disp_mat);
    ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);
    ElementData volumes        = m_model->compute_volumes();

    // Write results to the writer
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(disp_mat         , "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain           , "STRAIN");
    m_writer->write_eigen_matrix(stress           , "STRESS");
    m_writer->write_eigen_matrix(compliance_raw   , "COMPLIANCE_RAW");
    m_writer->write_eigen_matrix(compliance_adj   , "COMPLIANCE_ADJ");
    m_writer->write_eigen_matrix(dens_grad        , "DENS_GRAD");
    m_writer->write_eigen_matrix(volumes          , "VOLUME");

    // // build stiffness scalar
    // auto stiffness_scalar = density.array().pow(exponent);
    //
    // auto unconstrained = Timer::measure(
    //     [&]() { return this->m_model->build_unconstrained_index_matrix(); },
    //     "generating unconstrained index matrix"
    // );
    //
    // auto supp_vec = Timer::measure(
    //     [&]() { return this->m_model->build_support_vector(unconstrained, supps); },
    //     "building support vector"
    // );
    //
    // auto load_vec = Timer::measure(
    //     [&]() { return this->m_model->build_load_vector(unconstrained, loads); },
    //     "formulating load vector"
    // );
    //
    // auto stiffness = Timer::measure(
    //     [&]() { return this->m_model->build_stiffness_matrix(unconstrained, stiffness_scalar); },
    //     "constructing stiffness matrix"
    // );
    //
    // auto impl_load_vec = Timer::measure(
    //     [&]() { return this->m_model->build_implicit_load_vector(stiffness, supp_vec); },
    //     "computing implicit load vector"
    // );
    //
    // auto constrained = Timer::measure(
    //     [&]() { return this->m_model->build_constrained_index_matrix(supp_vec); },
    //     "creating constrained index matrix"
    // );
    //
    // auto mapping_vec = Timer::measure(
    //     [&]() { return this->m_model->build_mapping_vector(unconstrained, constrained); },
    //     "mapping unconstrained indices to constrained ones"
    // );
    //
    // load_vec += impl_load_vec;
    //
    // auto reduced_stiffness = Timer::measure(
    //     [&]() { return this->m_model->build_reduced_stiffness(mapping_vec, stiffness); },
    //     "reducing stiffness matrix"
    // );
    //
    // auto reduced_load = Timer::measure(
    //     [&]() { return this->m_model->build_reduced_load(mapping_vec, load_vec); },
    //     "reducing load vector"
    // );
    //
    // auto displacement = Timer::measure(
    //     [&]() { return solve(device, method, reduced_stiffness, reduced_load); },
    //     "solving system of equations"
    // );
    //
    // auto disp_matrix = m_model->build_global_displacement(constrained, displacement, unconstrained, supp_vec);
    //
    // NodeData stress;
    // NodeData strain;
    //
    // std::tie(stress, strain) = Timer::measure(
    //     [&]() { return m_model->compute_stress_strain(disp_matrix); },
    //     "Interpolation of stress and strain components to nodes"
    // );
    //
    // ElementData compliance_raw = m_model->compute_compliance(disp_matrix);
    // ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    // ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);
    // ElementData volumes        = m_model->compute_volumes();
    //
    // m_writer->add_loadcase(m_id);
    // m_writer->write_eigen_matrix(disp_matrix      , "DISPLACEMENT");
    // m_writer->write_eigen_matrix(strain           , "STRAIN");
    // m_writer->write_eigen_matrix(stress           , "STRESS");
    // m_writer->write_eigen_matrix(compliance_raw   , "COMPLIANCE_RAW");
    // m_writer->write_eigen_matrix(compliance_adj   , "COMPLIANCE_ADJ");
    // m_writer->write_eigen_matrix(dens_grad        , "DENS_GRAD");
    // m_writer->write_eigen_matrix(volumes          , "VOLUME");
}
