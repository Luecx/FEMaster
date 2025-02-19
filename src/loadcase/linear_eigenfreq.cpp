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


struct EigenMode {
    Precision eigenvalue;
    Precision eigenfreq;
    DynamicVector mode_shape;
    DynamicMatrix mode_shape_mat;
    Vec6 participation;

    EigenMode (Precision p_eigenvalue, DynamicVector p_mode_shape)
        : eigenvalue(p_eigenvalue),mode_shape(p_mode_shape) {
        eigenfreq = std::sqrt(eigenvalue) / (2 * M_PI);
    }

    void compute_participation(const DynamicMatrix& active_dof_vectors, const SparseMatrix& mass_matrix) {
        participation = Vec6::Zero();
        for (int i = 0; i < 6; i++) {
            participation(i) = active_dof_vectors.col(i).transpose() * (mass_matrix * mode_shape);
        }
    }

    void expand_mode_shape(int lagrange_dofs, const DynamicVector& active_lhs_vec, const IndexMatrix& active_dof_idx_mat) {
        DynamicVector mode_shape_active = mode_shape.head(mode_shape.size() - lagrange_dofs);

        auto expanded_mode_shape = mattools::expand_vec_to_vec(mode_shape_active, active_lhs_vec);

        mode_shape_mat = mattools::expand_vec_to_mat(active_dof_idx_mat, expanded_mode_shape);
    }

    // Comparator for sorting by eigenvalue
    bool operator<(const EigenMode &other) const {
        return eigenvalue < other.eigenvalue;
    }

};

// Helper function to pair eigenvalues and eigenvectors
std::vector<EigenMode> pair_eigenvalues_and_vectors(const DynamicVector &eigenvalues, const DynamicMatrix &mode_shapes) {
    std::vector<EigenMode> modes;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        modes.emplace_back(eigenvalues[i], mode_shapes.col(i));
    }
    return modes;
}

// Helper function to sort eigenvalue-vector pairs
void sort_eigen_pairs(std::vector<EigenMode> &modes) {
    std::sort(modes.begin(), modes.end());
}

int compute_scaling_exponent(const std::vector<EigenMode>& modes) {
    Precision max_coeff = 0;
    for (const EigenMode& mode : modes) {
        max_coeff = std::max(max_coeff, mode.participation.array().abs().maxCoeff());
    }
    Precision min_scaling = 1/max_coeff;
    Precision scaling_base = std::log10(min_scaling);

    // round it up
    scaling_base = std::ceil(scaling_base);

    return int(scaling_base);
}

// Function to display eigenvalues, eigenfrequencies, and participation to the console
void display_eigen_results(const std::vector<EigenMode> &modes) {
    int scaling_exponent = compute_scaling_exponent(modes);
    Precision s = std::pow(10, scaling_exponent);

    logging::info(true, std::setw(42), "", std::setw(33), "PARTICIPATION", " (x10^",scaling_exponent,")");
    logging::info(true,
                  std::setw(4), "Idx",
                  std::setw(20), "Eigenvalue",
                  std::setw(18), "Eigenfreq",
                  std::setw(12), "x", std::setw(8), "y", std::setw(8), "z",
                  std::setw(8), "rx", std::setw(8), "ry", std::setw(8), "rz");

    for (size_t i = 0; i < modes.size(); ++i) {
        const auto &mode = modes[i];
        logging::info(true, std::setw(4), i + 1,
                      std::setw(20), std::fixed, std::setprecision(6), Precision(mode.eigenvalue),
                      std::setw(18), std::fixed, std::setprecision(6), Precision(mode.eigenfreq),
                      std::setw(12), std::fixed, std::setprecision(3), Precision(mode.participation(0) * s),
                      std::setw(8) , std::fixed, std::setprecision(3), Precision(mode.participation(1) * s),
                      std::setw(8) , std::fixed, std::setprecision(3), Precision(mode.participation(2) * s),
                      std::setw(8) , std::fixed, std::setprecision(3), Precision(mode.participation(3) * s),
                      std::setw(8) , std::fixed, std::setprecision(3), Precision(mode.participation(4) * s),
                      std::setw(8) , std::fixed, std::setprecision(3), Precision(mode.participation(5) * s));
    }
}

void write_results(const std::vector<EigenMode>& modes, reader::Writer* writer, int loadcase_id) {
    writer->add_loadcase(loadcase_id);
    DynamicVector eigenvalues(modes.size());
    DynamicVector eigenfreqs(modes.size());
    DynamicMatrix mode_shapes(modes[0].mode_shape.size(), modes.size());
    DynamicMatrix participations(6, modes.size());

    for (size_t i = 0; i < modes.size(); ++i) {
        eigenvalues(i) = modes[i].eigenvalue;
        eigenfreqs(i) = modes[i].eigenfreq;

        writer->write_eigen_matrix(modes[i].mode_shape_mat, "MODE_SHAPE_" + std::to_string(i+1));
        writer->write_eigen_matrix(modes[i].participation, "PARTICIPATION_" + std::to_string(i+1));
    }

    writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    writer->write_eigen_matrix(DynamicMatrix(eigenfreqs), "EIGENFREQUENCIES");
    // writer->write_eigen_matrix(mode_shapes, "MODE_SHAPES");
    // writer->write_eigen_matrix(participations, "PARTICIPATIONS");
}


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
 * function to compute the vectors which have 1 at the active dofs for each respective axis
 * In total this is a Nx6 matrix where N is the number of active dofs. The first column is
 * 1 at the x-axis, the second at the y-axis and the third at the z-axis. The next three columns
 * are the same but for the rotation around the x, y and z axis.
 */
DynamicMatrix compute_active_dof_vectors(const IndexMatrix& active_dof_idx_mat, const int total_dofs, const DynamicVector& active_lhs_vec) {
    DynamicMatrix active_dof_vectors(total_dofs, 6);

    for (int i = 0; i < 6; i++) {
        DynamicMatrix sys_matrix = DynamicMatrix::Zero(active_dof_idx_mat.rows(), 6);
        // assign 1 to the respective column
        for(int j = 0; j < active_dof_idx_mat.rows(); j++) {
            sys_matrix(j, i) = 1;
        }

        // compute a vector
        auto active_dof_vector  = mattools::reduce_mat_to_vec(active_dof_idx_mat, sys_matrix);
        auto reduced_dof_vector = mattools::reduce_vec_to_vec(active_dof_vector, active_lhs_vec);

        // fill with zeros for lagrange multipliers
        int lagrange_dofs = total_dofs - reduced_dof_vector.size();
        reduced_dof_vector.conservativeResize(total_dofs);
        reduced_dof_vector.tail(lagrange_dofs) = DynamicVector::Zero(lagrange_dofs);

        // store in the active_dof_vectors matrix
        active_dof_vectors.col(i) = reduced_dof_vector;
    }

    return active_dof_vectors;
}


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

    m_model->assign_sections();

    // Step 1: Generate active_dof_idx_mat index matrix
    auto active_dof_idx_mat = Timer::measure(
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },  // Fixed method name
        "generating active_dof_idx_mat index matrix"
    );

    // Step 2: Build the global support matrix (includes all DOFs)
    auto global_supp_mat_pair = Timer::measure(
        [&]() { return this->m_model->build_support_matrix(supps); },  // Fixed parameter name
        "building global support matrix"
    );
    auto global_supp_mat = std::get<0>(global_supp_mat_pair);
    auto global_supp_eqs = std::get<1>(global_supp_mat_pair);

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
        [&]() { return m_model->build_constraint_matrix(active_dof_idx_mat, global_supp_eqs, characteristic_stiffness); },
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

    Index constrained_dofs   = active_dof_idx_mat.maxCoeff() + 1 + n - sol_rhs.rows();
    Index active_system_dofs = sol_rhs.rows() - n;

    Index n_eigenvalues = std::min((Index)num_eigenvalues, active_system_dofs);

    // Compress the stiffness matrix for efficient solving
    sol_mass_mat.makeCompressed();
    sol_stiffess_mat.makeCompressed();

    // Log system overview
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->_data->max_nodes);
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
                                               n_eigenvalues,
                                               true);
        },
        "solving constrained eigenvalue problem");

    // Step 8: Pair eigenvalues and eigenvectors, sort them, and compute participation
    auto eigenmodes = pair_eigenvalues_and_vectors(eigen_result.first, eigen_result.second);
    sort_eigen_pairs(eigenmodes);

    auto active_dof_vectors = compute_active_dof_vectors(active_dof_idx_mat, sol_rhs.rows(), active_lhs_vec);
    for(EigenMode& mode : eigenmodes) {
        mode.compute_participation(active_dof_vectors, sol_mass_mat);
        mode.expand_mode_shape(n, active_lhs_vec, active_dof_idx_mat);
    }
    display_eigen_results(eigenmodes);
    write_results(eigenmodes, m_writer, m_id);

    logging::info(true, "Eigenfrequency analysis completed successfully.");
}

} // namespace loadcase
} // namespace fem
