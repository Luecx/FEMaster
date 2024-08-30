/******************************************************************************
 * @file LinearStatic.cpp
 * @brief Implementation of the LinearStatic load case for performing linear
 * static analysis.
 *
 * @details This file implements the methods for performing linear static
 * analysis using the finite element method. The analysis computes displacements,
 * stresses, and reactions for a structure under static loading conditions. The
 * analysis assumes linear elasticity. The solver handles the formulation of
 * stiffness matrices, reduction of the system to account for boundary conditions,
 * and solution of the linear system of equations.
 *
 * The run method is responsible for assembling the system matrices and solving
 * for the displacements and other quantities of interest.
 *
 * @date Created on 04.09.2023
 ******************************************************************************/

#include "linear_static.h"

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_vec.h"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/extract_scaled_row_sum.h"

#include "../solve/eigen.h"

/**
 * @brief Constructs a LinearStatic load case with a given ID, writer, and model.
 *
 * @param id The unique ID of the load case.
 * @param writer Pointer to a Writer object for output handling.
 * @param model Pointer to the finite element model.
 */
fem::loadcase::LinearStatic::LinearStatic(ID id, reader::Writer* writer, model::Model* model)
    : LoadCase(id, writer, model) {}

/**
 * @brief Executes the linear static analysis, solving for displacements, stresses,
 * and reaction forces in the model.
 *
 * This method performs the following steps:
 * 1. Generate the unconstrained index matrix.
 * 2. Build the support vector based on boundary conditions.
 * 3. Formulate the load vector based on applied loads.
 * 4. Construct the stiffness matrix.
 * 5. Reduce the load vector and support vector to handle constraints.
 * 6. Solve the system of equations using the selected solver method.
 * 7. Expand the displacement vector to its full size and compute stresses and
 *    strains.
 * 8. Write the results (displacements, stresses, strains) to the output writer.
 */
void fem::loadcase::LinearStatic::run() {
    // Begin logging
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

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
        [&]() { return this->m_model->build_stiffness_matrix(unconstrained); },
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

    // Write results to the writer
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(disp_mat, "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain  , "STRAIN");
    m_writer->write_eigen_matrix(stress  , "STRESS");
}
