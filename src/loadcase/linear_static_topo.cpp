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
    : LinearStatic(id, writer, model), density(model->_data->max_elems, 1), orientation(model->_data->max_elems, 3) {
    density.setOnes();
    orientation.setZero();
}

void fem::loadcase::LinearStaticTopo::run() {

// Begin logging
    logging::info(true, "");
    logging::info(true, "");
    logging::info(true, "================================================================================================");
    logging::info(true, "LINEAR STATIC TOPO");
    logging::info(true, "================================================================================================");
    logging::info(true, "");

    auto stiffness_scalar = density.array().pow(exponent);

    m_model->_data->create_data(model::ElementDataEntries::TOPO_STIFFNESS, 1);
    m_model->_data->create_data(model::ElementDataEntries::TOPO_ANGLES   , 3);

    m_model->_data->get(model::ElementDataEntries::TOPO_STIFFNESS) = stiffness_scalar;
    m_model->_data->get(model::ElementDataEntries::TOPO_ANGLES)    = orientation;

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
        [&]() { return this->m_model->build_stiffness_matrix(active_dof_idx_mat, stiffness_scalar); },
        "constructing active stiffness matrix"
    );

    // compute characteristic stiffness by taking the mean of the diagonal
    Precision characteristic_stiffness = active_stiffness_mat.diagonal().mean();

    // Step 5: Construct the active Lagrangian constraint matrix
    auto active_lagrange_mat = Timer::measure(
        [&]() { return m_model->build_constraint_matrix(active_dof_idx_mat, characteristic_stiffness); },  // Fixed method name
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
    logging::info(true, "max nodes         : ", m_model->_data->max_nodes);
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

    ElementData compliance_raw = m_model->compute_compliance(global_disp_mat);
    ElementData compliance_adj = compliance_raw.array() * density.array().pow(exponent);
    ElementData dens_grad      = - exponent * compliance_raw.array() * density.array().pow(exponent - 1);
    ElementData volumes        = m_model->compute_volumes();
    ElementData angle_grad     = m_model->compute_compliance_angle_derivative(global_disp_mat);


    // Write results to the writer
    m_writer->add_loadcase(m_id);
    m_writer->write_eigen_matrix(global_disp_mat  , "DISPLACEMENT");
    m_writer->write_eigen_matrix(strain           , "STRAIN");
    m_writer->write_eigen_matrix(stress           , "STRESS");
    m_writer->write_eigen_matrix(compliance_raw   , "COMPLIANCE_RAW");
    m_writer->write_eigen_matrix(compliance_adj   , "COMPLIANCE_ADJ");
    m_writer->write_eigen_matrix(dens_grad        , "DENS_GRAD");
    m_writer->write_eigen_matrix(volumes          , "VOLUME");
    m_writer->write_eigen_matrix(density          , "DENSITY");
    m_writer->write_eigen_matrix(angle_grad       , "ORIENTATION_GRAD");
    m_writer->write_eigen_matrix(orientation      , "ORIENTATION");

    m_model->_data->remove(model::ElementDataEntries::TOPO_STIFFNESS);
    m_model->_data->remove(model::ElementDataEntries::TOPO_ANGLES);

}
