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
 * and reaction forces in the model using the finite element method (FEM).
 *
 * This method performs the following steps:
 * 1. **Generate the active DOF index matrix (`active_dof_idx_mat`)**:
 *    The degrees of freedom (DOFs) that are not constrained by boundary conditions
 *    are indexed, forming the `active_dof_idx_mat`. This matrix maps the active DOFs
 *    in the system, excluding those with prescribed displacements.
 *
 * 2. **Build the global support matrix for constraints (`global_supp_mat`)**:
 *    A matrix is built that represents the support conditions (boundary constraints)
 *    of the model. This matrix includes all DOFs (both active and constrained). The
 *    support matrix is used to ensure that the boundary conditions are applied correctly.
 *
 * 3. **Formulate the global load vector (`global_load_mat`)**:
 *    The load vector for the entire system is constructed based on the applied external
 *    forces. The load vector represents the forces acting on each degree of freedom
 *    (including both active and constrained DOFs).
 *
 * 4. **Construct the global stiffness matrix (`active_stiffness_mat`) and the Lagrangian constraint matrix (`active_lagrange_mat`)**:
 *    - The **stiffness matrix** represents the stiffness of the structure. It relates
 *      the displacements of the active DOFs to the forces acting on the structure.
 *      This matrix is generated for the active DOFs.
 *    - The **Lagrangian constraint matrix** (`C`), is formulated to enforce constraints
 *      using Lagrange multipliers. It contains the Lagrange multiplier equations that
 *      ensure the boundary conditions are satisfied.
 *
 * 5. **Assemble the full system matrix (`K C^T / C 0`)**:
 *    The global stiffness matrix and Lagrangian constraint matrix are combined into a
 *    larger system matrix. This matrix has the following block structure:
 *    \[
 *    \begin{bmatrix}
 *    K & C^T \\
 *    C & 0
 *    \end{bmatrix}
 *    \]
 *    Where:
 *    - `K` is the stiffness matrix for the active DOFs (size `m x m`, where `m` is the number of active DOFs).
 *    - `C^T` is the transpose of the Lagrangian constraint matrix (size `n x m`, where `n` is the number of Lagrangian multipliers).
 *    - `C` is the Lagrangian constraint matrix (size `n x m`).
 *    - The bottom-right block is a zero matrix of size `n x n` that corresponds to the Lagrangian multipliers.
 *
 * 6. **Reduce the load vector and support vector (`active_rhs_vec`, `active_lhs_vec`)**:
 *    The global load and support vectors are reduced to handle only the active DOFs.
 *    The active load vector (`active_rhs_vec`) contains the external forces acting
 *    on the active DOFs. Similarly, the active support vector (`active_lhs_vec`) contains
 *    the prescribed displacements for constrained DOFs (if any).
 *
 *    Additionally, the load vector is extended with Lagrangian multipliers to match
 *    the dimensions of the final system matrix.
 *
 * 7. **Solve the system of equations**:
 *    The final system of equations, represented by the full matrix
 *    \[
 *    \begin{bmatrix}
 *    K & C^T \\
 *    C & 0
 *    \end{bmatrix}
 *    \]
 *    is solved using the selected solver method (e.g., direct solver, iterative solver).
 *    The solution vector contains both the displacements of the active DOFs and the
 *    values of the Lagrangian multipliers that enforce the boundary conditions.
 *
 * 8. **Expand the displacement vector to its full size**:
 *    The reduced displacement vector (which includes only the active DOFs) is expanded
 *    back to the full size, including the constrained DOFs. This results in a displacement
 *    vector that corresponds to all the DOFs in the model.
 *
 * 9. **Compute stresses and strains**:
 *    Using the full displacement vector, stresses and strains are computed at the nodes
 *    of the model. These are derived from the displacement field using material properties
 *    and the finite element method (FEM).
 *
 * 10. **Write the results**:
 *    The computed displacements, stresses, and strains are written to the output writer,
 *    which stores the results for further processing or visualization. The output may
 *    include node-wise displacements, element-wise stresses and strains, and other relevant
 *    quantities.
 */
void fem::loadcase::LinearStatic::run() {
    // Begin logging
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC");
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
        [&]() { return this->m_model->build_load_matrix(loads); },
        "building global load matrix"
    );

    // Step 4: Construct the active stiffness matrix for active DOFs
    auto active_stiffness_mat = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing active stiffness matrix"
    );

    // Step 5: Construct the active Lagrangian constraint matrix
    auto active_lagrange_mat = Timer::measure(
        [&]() { return m_model->build_constraint_matrix(active_dof_idx_mat); },  // Fixed method name
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
            full_matrix.setFromTriplets(full_triplets.begin(), full_triplets.end());
            return full_matrix;
         },
        "assembling full lhs matrix including stiffness and Lagrangian"
    );

    auto active_lagrange_rhs = DynamicVector::Zero(n);  // Lagrangian RHS initialized to zero
    auto active_lagrange_lhs = DynamicVector::Constant(n, std::numeric_limits<Precision>::quiet_NaN());  // LHS for Lagrangian

    // Step 7: Reduce global load and support matrices to active DOFs
    auto active_rhs_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat); },
        "reducing load vector (RHS)"
    );
    auto active_lhs_vec = Timer::measure(
        [&]() { return mattools::reduce_mat_to_vec(active_dof_idx_mat, global_supp_mat); },
        "reducing support vector (LHS)"
    );

    // Extend RHS and LHS vectors for Lagrangian DOFs
    DynamicVector full_rhs_vec(m + n);
    DynamicVector full_lhs_vec(m + n);
    full_rhs_vec << active_rhs_vec, active_lagrange_rhs;  // Combine active RHS with Lagrangian RHS
    full_lhs_vec << active_lhs_vec, active_lagrange_lhs;  // Combine active LHS with Lagrangian LHS

    // Step 8: Compute the implicit load vector to account for constraints
    auto implicit_rhs_vec = Timer::measure(
        [&]() { return mattools::extract_scaled_row_sum(active_lhs_mat, full_lhs_vec); },
        "computing implicit load vector"
    );

    // Add implicit load vector to the RHS
    full_rhs_vec += implicit_rhs_vec;

    // Step 9: Reduce full system matrix and RHS vector to handle constrained DOFs
    auto sol_matrix = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(active_lhs_mat, full_lhs_vec); },
        "reducing stiffness matrix to solver-ready form"
    );
    auto sol_rhs = Timer::measure(
        [&]() { return mattools::reduce_vec_to_vec(full_rhs_vec, full_lhs_vec); },
        "reducing load vector (RHS) to solver-ready form"
    );

    // Compress the stiffness matrix for efficient solving
    sol_matrix.makeCompressed();

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

    // Step 10: Solve the system of equations to get displacements
    auto sol_lhs = Timer::measure(
        [&]() { return solve(device, method, sol_matrix, sol_rhs); },
        "solving system of equations"
    );

    // Step 11: Expand the reduced displacement vector to full size
    auto active_disp_vec = Timer::measure(
        [&]() {
            DynamicVector active_disp_u = sol_lhs.segment(0, sol_lhs.rows() - n);  // Extract active DOF displacements
            return mattools::expand_vec_to_vec(active_disp_u, active_lhs_vec);     // Expand back to full displacement
        },
        "expanding displacement vector"
    );

    // Expand the displacement vector to global matrix form
    auto global_disp_mat = Timer::measure(
        [&]() { return mattools::expand_vec_to_mat(active_dof_idx_mat, active_disp_vec); },
        "expanding displacement vector to matrix form"
    );

    // Step 12: Compute stresses and strains at the nodes
    NodeData stress;
    NodeData strain;

    std::tie(stress, strain) = Timer::measure(
        [&]() { return m_model->compute_stress_strain(global_disp_mat); },
        "Interpolating stress and strain at nodes"
    );

    // Write results to the writer
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(global_disp_mat, "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain, "STRAIN");
    m_writer->write_eigen_matrix(stress, "STRESS");
}
