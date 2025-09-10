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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>

#include "../mattools/assemble.tpp"
#include "../mattools/reduce_mat_to_mat.h"
#include "../mattools/reduce_vec_to_vec.h"
#include "../mattools/reduce_mat_to_vec.h"

#include "../solve/eigen.h"  // uses new API: eigs(...) returning vector<EigenValueVectorPair>

namespace fem {
namespace loadcase {

// Result container for one eigenmode (value, freq, vectors, etc.)
struct EigenMode {
    Precision     eigenvalue;    // λ
    Precision     eigenfreq;     // f = sqrt(λ)/(2π)
    DynamicVector mode_shape;    // reduced (contains Lagrange rows at tail)
    DynamicMatrix mode_shape_mat;// expanded to node x dof layout
    Vec6          participation; // modal participation (x,y,z,rx,ry,rz)

    EigenMode(Precision p_eigenvalue, DynamicVector p_mode_shape)
        : eigenvalue(p_eigenvalue), mode_shape(std::move(p_mode_shape))
    {
        eigenfreq = std::sqrt(eigenvalue) / (2 * M_PI);
    }

    void compute_participation(const DynamicMatrix& active_dof_vectors,
                               const SparseMatrix&  mass_matrix)
    {
        participation = Vec6::Zero();
        for (int i = 0; i < 6; i++) {
            participation(i) = active_dof_vectors.col(i).transpose() * (mass_matrix * mode_shape);
        }
    }

    void expand_mode_shape(int lagrange_dofs,
                           const DynamicVector& active_lhs_vec,
                           const IndexMatrix&   active_dof_idx_mat)
    {
        DynamicVector mode_shape_active = mode_shape.head(mode_shape.size() - lagrange_dofs);
        auto expanded_mode_shape        = mattools::expand_vec_to_vec(mode_shape_active, active_lhs_vec);
        mode_shape_mat                  = mattools::expand_vec_to_mat(active_dof_idx_mat, expanded_mode_shape);
    }

    bool operator<(const EigenMode& other) const {
        return eigenvalue < other.eigenvalue;
    }
};

// Helper: convert returned eigenpairs to EigenMode list
static std::vector<EigenMode>
    make_modes_from_pairs(const std::vector<solver::EigenValueVectorPair>& pairs)
{
    std::vector<EigenMode> modes;
    modes.reserve(pairs.size());
    for (const auto& p : pairs) {
        modes.emplace_back(p.value, p.vector);
    }
    return modes;
}

// Helper: sort modes by eigenvalue
static void sort_eigen_pairs(std::vector<EigenMode>& modes) {
    std::sort(modes.begin(), modes.end());
}

static int compute_scaling_exponent(const std::vector<EigenMode>& modes) {
    Precision max_coeff = 0;
    for (const EigenMode& mode : modes) {
        max_coeff = std::max(max_coeff, mode.participation.array().abs().maxCoeff());
    }
    if (max_coeff <= Precision(0)) return 0;
    Precision min_scaling  = Precision(1) / max_coeff;
    Precision scaling_base = std::log10(min_scaling);
    return static_cast<int>(std::ceil(scaling_base));
}

// Function to display eigenvalues, eigenfrequencies, and participation to the console
static void display_eigen_results(const std::vector<EigenMode>& modes) {
    const int scaling_exponent = compute_scaling_exponent(modes);
    const Precision s = std::pow(Precision(10), scaling_exponent);

    logging::info(true, std::setw(42), "", std::setw(33), "PARTICIPATION", " (x10^", scaling_exponent, ")");
    logging::info(true,
                  std::setw(4),  "Idx",
                  std::setw(20), "Eigenvalue",
                  std::setw(18), "Eigenfreq",
                  std::setw(12), "x",  std::setw(8), "y",  std::setw(8), "z",
                  std::setw(8),  "rx", std::setw(8), "ry", std::setw(8), "rz");

    for (size_t i = 0; i < modes.size(); ++i) {
        const auto& mode = modes[i];
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

static void write_results(const std::vector<EigenMode>& modes,
                          reader::Writer*               writer,
                          int                           loadcase_id)
{
    writer->add_loadcase(loadcase_id);

    DynamicVector eigenvalues(modes.size());
    DynamicVector eigenfreqs (modes.size());

    for (size_t i = 0; i < modes.size(); ++i) {
        eigenvalues(i) = modes[i].eigenvalue;
        eigenfreqs (i) = modes[i].eigenfreq;

        writer->write_eigen_matrix(modes[i].mode_shape_mat, "MODE_SHAPE_"   + std::to_string(i + 1));
        writer->write_eigen_matrix(modes[i].participation,  "PARTICIPATION_" + std::to_string(i + 1));
    }

    writer->write_eigen_matrix(DynamicMatrix(eigenvalues), "EIGENVALUES");
    writer->write_eigen_matrix(DynamicMatrix(eigenfreqs ), "EIGENFREQUENCIES");
    // If you want stacked matrices of all modes/participations, build and write them here.
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
fem::loadcase::LinearEigenfrequency::LinearEigenfrequency(ID id,
                                                          reader::Writer* writer,
                                                          model::Model*   model,
                                                          int             numEigenvalues)
    : LoadCase(id, writer, model), num_eigenvalues(numEigenvalues) {}

/**
 * @brief Computes the axis unit vectors restricted to active DOFs (and extended with Lagrange rows).
 *
 * Constructs a matrix of size (total_dofs x 6) where columns are unit load-like
 * vectors in the translational (x,y,z) and rotational (rx,ry,rz) directions for
 * the current active system, including zero-padding for Lagrange multipliers.
 */
static DynamicMatrix
    compute_active_dof_vectors(const IndexMatrix&  active_dof_idx_mat,
                               const int           total_dofs,
                               const DynamicVector& active_lhs_vec)
{
    DynamicMatrix active_dof_vectors(total_dofs, 6);

    for (int i = 0; i < 6; i++) {
        DynamicMatrix sys_matrix = DynamicMatrix::Zero(active_dof_idx_mat.rows(), 6);
        for (int j = 0; j < active_dof_idx_mat.rows(); j++) {
            sys_matrix(j, i) = 1;
        }

        auto active_dof_vector  = mattools::reduce_mat_to_vec(active_dof_idx_mat, sys_matrix);
        auto reduced_dof_vector = mattools::reduce_vec_to_vec(active_dof_vector, active_lhs_vec);

        // Fill zeros for Lagrange multipliers
        const int lagrange_dofs = total_dofs - static_cast<int>(reduced_dof_vector.size());
        reduced_dof_vector.conservativeResize(total_dofs);
        reduced_dof_vector.tail(lagrange_dofs) = DynamicVector::Zero(lagrange_dofs);

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
        [&]() { return this->m_model->build_unconstrained_index_matrix(); },
        "generating active_dof_idx_mat index matrix"
    );

    // Step 2: Build the global support matrix (includes all DOFs)
    auto global_supp_mat_pair = Timer::measure(
        [&]() { return this->m_model->build_support_matrix(supps); },
        "building global support matrix"
    );
    auto global_supp_mat = std::get<0>(global_supp_mat_pair);
    auto global_supp_eqs = std::get<1>(global_supp_mat_pair);

    // Step 3: Build the global load matrix based on applied loads (includes all DOFs)
    auto global_load_mat = Timer::measure(
        [&]() { return this->m_model->build_load_matrix(); },
        "building global load matrix"
    );

    // Step 4: Construct the active stiffness and mass matrices for active DOFs
    auto active_stiffness_mat = Timer::measure(
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat); },
        "constructing active stiffness matrix"
    );
    auto active_mass_mat = Timer::measure(
        [&]() { return this->m_model->build_lumped_mass_matrix(active_dof_idx_mat); },
        "constructing mass matrix"
    );

    // Characteristic stiffness (for Lagrange regularization)
    const Precision characteristic_stiffness = active_stiffness_mat.diagonal().mean();

    // Step 5: Construct the active Lagrangian constraint matrix
    auto active_lagrange_mat = Timer::measure(
        [&]() { return m_model->build_constraint_matrix(active_dof_idx_mat,
                                                        global_supp_eqs,
                                                        characteristic_stiffness); },
        "constructing active Lagrangian matrix"
    );

    const int m   = static_cast<int>(active_stiffness_mat.rows()); // active DOFs
    const int n   = static_cast<int>(active_lagrange_mat.rows());  // # Lagrange multipliers
    const int nnz = static_cast<int>(active_stiffness_mat.nonZeros()
                                     + 2 * active_lagrange_mat.nonZeros());

    // Step 6: Assemble the full system matrix (stiffness + Lagrangian)
    auto active_lhs_mat = Timer::measure(
        [&]() {
            SparseMatrix full_matrix(m + n, m + n);
            TripletList full_triplets;
            full_triplets.reserve(nnz);

            // Insert stiffness matrix
            for (int k = 0; k < active_stiffness_mat.outerSize(); ++k) {
                for (SparseMatrix::InnerIterator it(active_stiffness_mat, k); it; ++it) {
                    full_triplets.emplace_back(it.row(), it.col(), it.value());
                }
            }

            // Insert Lagrange blocks
            for (int k = 0; k < active_lagrange_mat.outerSize(); ++k) {
                for (SparseMatrix::InnerIterator it(active_lagrange_mat, k); it; ++it) {
                    full_triplets.emplace_back(it.row() + m, it.col(), it.value());
                    full_triplets.emplace_back(it.col(), it.row() + m, it.value());
                }
            }

            // Regularization (small negative on Λ-Λ block)
            for (int i = 0; i < n; i++) {
                full_triplets.emplace_back(m + i, m + i, - characteristic_stiffness / 1e6);
            }

            full_matrix.setFromTriplets(full_triplets.begin(), full_triplets.end());
            return full_matrix;
        },
        "assembling full lhs matrix including stiffness and Lagrangian"
    );

    // Resize the mass matrix to (m+n) to align with augmented system
    active_mass_mat.conservativeResize(m + n, m + n);

    // Lagrange RHS/LHS
    auto active_lagrange_rhs = DynamicVector::Zero(n);
    auto active_lagrange_lhs = DynamicVector::Constant(n, std::numeric_limits<Precision>::quiet_NaN());

    // Reduce global load/support to active vectors
    auto active_rhs_vec = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_load_mat);
    auto active_lhs_vec = mattools::reduce_mat_to_vec(active_dof_idx_mat, global_supp_mat);

    // Extend RHS/LHS by Lagrange DOFs
    DynamicVector full_rhs_vec(m + n);
    DynamicVector full_lhs_vec(m + n);
    full_rhs_vec << active_rhs_vec, active_lagrange_rhs;
    full_lhs_vec << active_lhs_vec, active_lagrange_lhs;

    auto sol_rhs = mattools::reduce_vec_to_vec(full_rhs_vec, full_lhs_vec);

    // Step 9: Reduce augmented matrices for constrained DOFs
    auto sol_stiffess_mat = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(active_lhs_mat, full_lhs_vec); },
        "reducing stiffness matrix to solver-ready form"
    );
    auto sol_mass_mat = Timer::measure(
        [&]() { return mattools::reduce_mat_to_mat(active_mass_mat, full_lhs_vec); },
        "reducing mass matrix to solver-ready form"
    );

    const Index constrained_dofs   = active_dof_idx_mat.maxCoeff() + 1 + n - sol_rhs.rows();
    const Index active_system_dofs = sol_rhs.rows() - n;

    Index n_eigenvalues = std::min<Index>(static_cast<Index>(num_eigenvalues), active_system_dofs);

    // Compress for Spectra/Eigen
    sol_mass_mat.makeCompressed();
    sol_stiffess_mat.makeCompressed();

    // Overview
    logging::info(true, "");
    logging::info(true, "Overview");
    logging::up();
    logging::info(true, "max nodes         : ", m_model->_data->max_nodes);
    logging::info(true, "system total DOFs : ", active_dof_idx_mat.maxCoeff() + 1);
    logging::info(true, "lagrange DOFs     : ", n);
    logging::info(true, "total DOFs        : ", active_dof_idx_mat.maxCoeff() + 1 + n);
    logging::info(true, "constrained DOFs  : ", constrained_dofs);
    logging::info(true, "final DOFs        : ", sol_rhs.rows());
    logging::down();

    // Step 7: Solve constrained generalized eigenvalue problem K x = λ M x
    // Use new eigen API: always returns pairs (value + vector). For eigenfrequencies,
    // choose ShiftInvert with σ=0 and sorting SmallestAlge to get the lowest ω^2 first.
    auto eigen_pairs = Timer::measure(
        [&]() {
            solver::EigenOpts opts;
            opts.mode = solver::EigenMode::ShiftInvert;          // mass PD → shift-invert
            opts.sigma = 0.0;                                    // around zero
            opts.sort = solver::EigenOpts::Sort::LargestMagn;    // smallest λ = ω^2 first
            return solver::eigs(solver::CPU, sol_stiffess_mat, sol_mass_mat,
                                static_cast<int>(n_eigenvalues), opts);
        },
        "solving constrained eigenvalue problem"
    );

    // Step 8: Build EigenMode objects, sort, compute participation, expand shapes
    auto eigenmodes = make_modes_from_pairs(eigen_pairs);
    sort_eigen_pairs(eigenmodes);

    auto active_dof_vectors = compute_active_dof_vectors(active_dof_idx_mat, static_cast<int>(sol_rhs.rows()), active_lhs_vec);
    for (EigenMode& mode : eigenmodes) {
        mode.compute_participation(active_dof_vectors, sol_mass_mat);
        mode.expand_mode_shape(n, active_lhs_vec, active_dof_idx_mat);
    }

    display_eigen_results(eigenmodes);
    write_results(eigenmodes, m_writer, m_id);

    logging::info(true, "Eigenfrequency analysis completed successfully.");
}

} // namespace loadcase
} // namespace fem
