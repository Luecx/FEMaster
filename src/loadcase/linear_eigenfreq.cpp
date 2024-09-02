/******************************************************************************
 * @file LinearEigenfrequency.cpp
 * @brief Implementation of the LinearEigenfrequency load case for performing
 * linear eigenfrequency analysis.
 *
 * @details This file implements the methods for performing linear eigenfrequency
 * analysis using the finite element method. The analysis computes the natural
 * frequencies and mode shapes of a structure under free vibration conditions.
 * It solves the generalized eigenvalue problem involving the stiffness and mass
 * matrices.
 *
 * The run method assembles the stiffness and mass matrices, reduces them to account
 * for boundary conditions, and solves the eigenvalue problem.
 *
 * @date Created on <Your Creation Date>
 ******************************************************************************/

#include "linear_eigenfreq.h"

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/reduce_mat_to_vec.h"

#include "../solve/eigen.h"

namespace fem {
namespace loadcase {

/**
 * @brief Constructs a LinearEigenfrequency load case with a given ID, writer, model,
 * and the number of eigenvalues to compute.
 *
 * @param id The unique ID of the load case.
 * @param writer Pointer to a Writer object for output handling.
 * @param model Pointer to the finite element model.
 * @param numEigenvalues The number of eigenvalues (natural frequencies) to compute.
 */
fem::loadcase::LinearEigenfrequency::LinearEigenfrequency(ID id, reader::Writer* writer, model::Model* model, int numEigenvalues)
    : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

/**
 * @brief Executes the linear eigenfrequency analysis, solving for natural frequencies
 * and mode shapes of the model.
 *
 * This method performs the following steps:
 * 1. Generate the unconstrained index matrix.
 * 2. Build the support vector based on boundary conditions.
 * 3. Construct the stiffness matrix and mass matrix.
 * 4. Reduce the matrices to handle constraints.
 * 5. Solve the eigenvalue problem for the desired number of eigenvalues.
 * 6. Write the mode shapes to the output writer.
 */
void fem::loadcase::LinearEigenfrequency::run() {
    // Begin logging
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR EIGENFREQUENCY");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    // Step 1: Generate unconstrained index matrix
    auto unconstrained = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating unconstrained index matrix"
    );

    // Step 2: Build the support matrix for the constraints
    auto supp_mat = Timer::measure(
        [&]() { return this->m_model->build_constraint_matrix(supps); },
        "building support vector"
    );

    // Step 3: Construct the global stiffness matrix
    auto stiffness = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(unconstrained); },
        "constructing stiffness matrix"
    );

    // Step 3b: Construct the global lumped mass matrix
    auto mass = Timer::measure(
        [&]() { return this->m_model->build_lumped_mass_matrix(unconstrained); },
        "constructing lumped mass matrix"
    );

    // Step 4: Reduce stiffness matrix and mass matrix to handle constrained degrees of freedom
    auto reduced_supp_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(unconstrained, supp_mat); },
        "reducing support vector"
    );

    auto sol_stiffness = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(stiffness, reduced_supp_vec); },
        "reducing stiffness matrix"
    );

    auto sol_mass = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(mass, reduced_supp_vec); },
        "reducing mass matrix"
    );
    sol_stiffness.makeCompressed();
    sol_mass.makeCompressed();


    // Log the overview of the system
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->max_nodes);
    logging::info(true, "total dofs        : ", unconstrained.maxCoeff() + 1);
    logging::info(true, "unconstrained dofs: ", sol_stiffness.rows());
    logging::info(true, "constrained dofs  : ", unconstrained.maxCoeff() + 1 - sol_stiffness.rows());
    logging::down();

    // Step 5: Solve the eigenvalue problem
    auto eigen_result = Timer::measure(
        [&]() { return solver::compute_eigenvalues(solver::CPU, sol_stiffness, sol_mass, num_eigenvalues, true); },
        "solving eigenvalue problem"
    );

    DynamicVector eigenvalues = eigen_result.first;
    DynamicVector eigenfreqs  = eigenvalues.array().abs().sqrt() / (2 * M_PI);
    DynamicMatrix mode_shapes = eigen_result.second;

    // Log the eigenvalues and eigenfrequencies
    for (int i = 0; i < num_eigenvalues; i++) {
        logging::info(true, "Eigenvalue ", std::setw(4), i+1, " : ",
                      std::setw(20), std::fixed, Precision(eigenvalues(i)), " ,",
                      std::setw(20), std::fixed, Precision(eigenfreqs(i)), " cycles / time period");
    }

    // write results
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    m_writer->write_eigen_matrix(DynamicMatrix(eigenfreqs), "EIGENFREQUENCIES");
    for (int i = 0; i < num_eigenvalues; i++) {
        DynamicVector mode_shape = mode_shapes.col(i);
        // extend mode shape to full size
        // Step 8: Expand the reduced displacement vector to full size
        auto shape_red = mattools::expand_vec_to_vec(mode_shape, reduced_supp_vec);
        auto shape_mat = mattools::expand_vec_to_mat(unconstrained, shape_red);

        m_writer->write_eigen_matrix(shape_mat, "MODE_SHAPE_" + std::to_string(i));
    }
    // Step 6: Write the mode shapes to the writer
    //m_writer->add_loadcase(m_id);
    //m_writer->write_eigen_matrix(mode_shapes, "MODE_SHAPES");

    // Log the completion of the eigenfrequency analysis
    logging::info(true, "Eigenfrequency analysis completed successfully.");
}

} // namespace loadcase
} // namespace fem
