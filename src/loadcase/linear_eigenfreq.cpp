/******************************************************************************
* @file LinearEigenfrequency.cpp
* @brief Implementation of the LinearEigenfrequency load case for performing
* linear eigenfrequency analysis.
*
* @details This file implements the methods for performing linear eigenfrequency
* analysis using the finite element method. The analysis computes the natural
* frequencies and mode shapes of a structure under free vibration conditions.
* It solves the generalized eigenvalue problem involving the stiffness and mass
* matrices with Lagrange multiplier constraints.
*
* The run method assembles the stiffness, mass, and constraint matrices, reduces them
* to account for boundary conditions, and solves the constrained eigenvalue problem.
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
* 1. Generate the active DOF index matrix.
* 2. Build the support matrix based on boundary conditions.
* 3. Construct the stiffness matrix, mass matrix, and Lagrangian matrix.
* 4. Assemble the augmented system matrices.
* 5. Apply Lagrangian multipliers for constraints.
* 6. Solve the constrained eigenvalue problem.
* 7. Write the mode shapes to the output writer.
 */
void fem::loadcase::LinearEigenfrequency::run() {

    // Begin logging
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR EIGENFREQUENCY");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    // Step 1: Generate active_dof_idx_mat index matrix
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },  // Fixed method name
        "generating active_dof_idx_mat index matrix"
    );

    // Step 2: Build the global support matrix (includes all DOFs)
    auto global_supp_mat = Timer::measure(
        [&]() { return this->m_model->build_support_matrix(supps); },  // Fixed parameter name
        "building global support matrix"
    );

    // Step 3: Build the global load matrix based on applied loads (includes all DOFs)
    auto global_load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(); },
        "building global load matrix"
    );

    // Step 4: Construct the active stiffness matrix for active DOFs
    auto active_stiffness_mat = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing active stiffness matrix"
    );
    auto active_mass_mat = Timer::measure(
        [&]() { return this->m_model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix"
    );

    // compute characteristic stiffness by taking the mean of the diagonal
    Precision characteristic_stiffness = active_stiffness_mat.diagonal().mean();

    // Step 5: Construct the active Lagrangian constraint matrix
    auto active_lagrange_mat = Timer::measure(
        [&]() { return m_model->build_constraint_matrix(active_dof_idx_mat, characteristic_stiffness); },
        "constructing active Lagrangian matrix"
    );

    int m   = active_stiffness_mat.rows();                   // Number of active DOFs
    int n   = active_lagrange_mat.rows();                    // Number of Lagrangian multipliers
    int nnz = active_stiffness_mat.nonZeros() + 2 * active_lagrange_mat.nonZeros();

    // Step 6: Assemble the full system matrix (stiffness + Lagrangian)
    auto active_lhs_mat = Timer::measure(
        [&]() {
            SparseMatrix full_matrix(m + n, m + n);
            TripletList full_triplets;
            full_triplets.reserve(nnz);

            // Insert stiffness matrix into full system matrix
            for (int k = 0; k < active_stiffness_mat.outerSize(); ++k) {
                for (SparseMatrix::InnerIterator it(active_stiffness_mat, k); it; ++it) {
                    full_triplets.push_back(Triplet(it.row(), it.col(), it.value()));
                }
            }

            // Insert Lagrangian matrix into full system matrix
            for (int k = 0; k < active_lagrange_mat.outerSize(); ++k) {
                for (SparseMatrix::InnerIterator it(active_lagrange_mat, k); it; ++it) {
                    full_triplets.push_back(Triplet(it.row() + m, it.col(), it.value()));
                    full_triplets.push_back(Triplet(it.col(), it.row() + m, it.value()));
                }
            }

            // insert regularization term at the bottom right
            for (int i = 0; i < n; i++) {
                full_triplets.push_back(Triplet(m + i, m + i, - characteristic_stiffness / 1e6));
            }
            full_matrix.setFromTriplets(full_triplets.begin(), full_triplets.end());
            return full_matrix;
        },
        "assembling full lhs matrix including stiffness and Lagrangian"
    );
    // resize the mass matrix
    active_mass_mat.conservativeResize(m + n, m + n);

    auto active_lagrange_rhs = DynamicVector::Zero(n);  // Lagrangian RHS initialized to zero
    auto active_lagrange_lhs = DynamicVector::Constant(n, std::numeric_limits<Precision>::quiet_NaN());  // LHS for Lagrangian

    auto active_rhs_vec = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat);
    auto active_lhs_vec = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_supp_mat);

    // Extend RHS and LHS vectors for Lagrangian DOFs
    DynamicVector full_rhs_vec(m + n);
    DynamicVector full_lhs_vec(m + n);
    full_rhs_vec << active_rhs_vec, active_lagrange_rhs;  // Combine active RHS with Lagrangian RHS
    full_lhs_vec << active_lhs_vec, active_lagrange_lhs;  // Combine active LHS with Lagrangian LHS

    auto sol_rhs = mattools::reduce_vec_to_vec(full_rhs_vec, full_lhs_vec);

    // Step 9: Reduce full system matrix and RHS vector to handle constrained DOFs
    auto sol_stiffess_mat = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(active_lhs_mat, full_lhs_vec); },
        "reducing stiffness matrix to solver-ready form"
    );
    auto sol_mass_mat = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(active_mass_mat, full_lhs_vec); },
        "reducing mass matrix to solver-ready form"
    );


    // Compress the stiffness matrix for efficient solving
    sol_mass_mat.makeCompressed();
    sol_stiffess_mat.makeCompressed();

    // Log system overview
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->max_nodes);
    logging::info(true, "system total DOFs : ", active_dof_idx_mat.maxCoeff() + 1);
    logging::info(true, "lagrange DOFs     : ", n);
    logging::info(true, "total DOFs        : ", active_dof_idx_mat.maxCoeff() + 1 + n);
    logging::info(true, "constrained DOFs  : ", active_dof_idx_mat.maxCoeff() + 1 + n - sol_rhs.rows());
    logging::info(true, "final DOFs        : ", sol_rhs.rows());
    logging::down();


    // Step 7: Solve the constrained eigenvalue problem
    auto eigen_result = Timer::measure(
        [&]() {
            return solver::compute_eigenvalues(solver::CPU,
                                               sol_stiffess_mat,
                                               sol_mass_mat,
                                               num_eigenvalues,
                                               true);
        },
        "solving constrained eigenvalue problem");

    DynamicVector eigenvalues = eigen_result.first;
    DynamicVector eigenfreqs = eigenvalues.array().abs().sqrt() / (2 * M_PI);
    DynamicMatrix mode_shapes = eigen_result.second;

    // Log the eigenvalues and eigenfrequencies
    for (int i = 0; i < num_eigenvalues; i++) {
        logging::info(true, "Eigenvalue ", std::setw(4), i + 1, " : ", std::setw(20), std::fixed, Precision(eigenvalues(i)), " ,", std::setw(20), std::fixed, Precision(eigenfreqs(i)), " cycles / time period");
    }

    // Write results
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    m_writer->write_eigen_matrix(DynamicMatrix(eigenfreqs), "EIGENFREQUENCIES");
    for (int i = 0; i < num_eigenvalues; i++) {
        // Extract mode shape without Lagrange multipliers
        DynamicVector mode_shape_active = mode_shapes.col(i).head(m);

        // Expand mode shape to full displacement vector size
        auto expanded_mode_shape = mattools::expand_vec_to_vec(mode_shape_active, active_lhs_vec);

        // Convert to matrix form for output
        auto shape_mat = mattools::expand_vec_to_mat(active_dof_idx_mat, expanded_mode_shape);

        // Write the expanded mode shape to file
        m_writer->write_eigen_matrix(shape_mat, "MODE_SHAPE_" + std::to_string(i));
    }

    // Log the completion of the eigenfrequency analysis
    logging::info(true, "Eigenfrequency analysis completed successfully.");
}

} // namespace loadcase
} // namespace fem
