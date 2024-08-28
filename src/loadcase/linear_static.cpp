//
// Created by Luecx on 04.09.2023.
//
#include "linear_static.h"

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"

#include "../solve/eigen.h"

fem::loadcase::LinearStatic::LinearStatic(ID id, reader::Writer* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

void fem::loadcase::LinearStatic::run() {
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    auto unconstrained = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating unconstrained index matrix"
    );

    auto supp_mat = Timer::measure(
        [&]() { return this->m_model->build_constraint_matrix(supps);},
        "building support vector"
    );

    auto load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(loads); },
        "formulating load vector"
    );

    auto stiffness = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(unconstrained); },
        "constructing stiffness matrix"
    );

    auto reduced_load_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(unconstrained, load_mat); },
        "reducing load vector"
    );

    auto reduced_supp_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(unconstrained, supp_mat); },
        "reducing support vector"
    );

    auto impl_load_vec = Timer::measure(
        [&]() { return mattools::extract_scaled_row_sum(stiffness, reduced_supp_vec); },
        "computing implicit load vector"
    );

    reduced_load_vec += impl_load_vec;

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
    sol_stiffness.makeCompressed();

    // show overview
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->max_nodes);
    logging::info(true, "total dofs        : ", unconstrained.maxCoeff() + 1);
    logging::info(true, "unconstrained dofs: ", sol_support.rows());
    logging::info(true, "constrained dofs  : ", unconstrained.maxCoeff() + 1 - sol_support.rows());
    logging::down();

    auto sol_disp = Timer::measure(
        [&]() { return solve(device, method, sol_stiffness, sol_load); },
        "solving system of equations"
    );

    // expand to reduced size first
    auto reduced_disp = Timer::measure(
        [&]() { return mattools::expand_vec_to_vec(sol_disp, reduced_supp_vec); },
        "expanding displacement vector"
    );

    // expand to full size
    auto disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(unconstrained, reduced_disp); },
        "expanding displacement vector"
    );

    NodeData stress;
    NodeData strain;

    std::tie(stress, strain) = Timer::measure(
        [&]() { return m_model->compute_stress_strain(disp_mat); },
        "Interpolation of stress and strain components to nodes"
    );

    // sol_stiffness.makeCompressed();
    // std::pair<DynamicVector, DynamicMatrix> res = solver::compute_eigenvalues(solver::CPU, sol_stiffness, 10, false);
    // std::cout << res.first << std::endl;

    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(disp_mat, "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain  , "STRAIN");
    m_writer->write_eigen_matrix(stress  , "STRESS");
}
