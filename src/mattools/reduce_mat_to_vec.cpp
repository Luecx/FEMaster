/******************************************************************************
 * @file reduce_mat_to_vec.cpp
 * @brief Implements functions to reduce a SystemDofIds matrix and NodeData into
 * a DynamicVector and to expand it back.
 *
 * The reduction ignores inactive DOFs marked by -1 in the SystemDofIds, while
 * the expansion fills the NodeData based on the active DOFs and the given reduced
 * DynamicVector.
 *
 * @date Created on 28.08.2024
 ******************************************************************************/

#include "reduce_mat_to_vec.h"
#include <iostream>  // Optional: For debug output
#include "../core/core.h"

namespace fem { namespace mattools {

DynamicVector reduce_mat_to_vec(const SystemDofIds& dof_ids, const NodeData& bc_matrix) {
    std::vector<Precision> active_dofs;

    // Iterate through all rows and columns of the SystemDofIds matrix
    for (int i = 0; i < dof_ids.rows(); ++i) {
        for (int j = 0; j < dof_ids.cols(); ++j) {
            // Check if the DOF is active (i.e., not -1)
            if (dof_ids(i, j) != -1) {
                active_dofs.push_back(bc_matrix(i, j));
            }
        }
    }

    // Copy active DOFs into a DynamicVector (Eigen::VectorXd)
    DynamicVector reduced_vector(active_dofs.size());
    for (Index i = 0; i < active_dofs.size(); ++i) {
        reduced_vector(i) = active_dofs[i];
    }

    return reduced_vector;
}

NodeData expand_vec_to_mat(const SystemDofIds& dof_ids, const DynamicVector& reduced_vector) {
    NodeData expanded_matrix = NodeData::Zero(dof_ids.rows(), dof_ids.cols());
    int reduced_index = 0;

    // Iterate through all rows and columns of the SystemDofIds matrix
    for (int i = 0; i < dof_ids.rows(); ++i) {
        for (int j = 0; j < dof_ids.cols(); ++j) {
            // If the DOF is active (not -1), assign the value from the reduced vector
            if (dof_ids(i, j) != -1) {
                expanded_matrix(i, j) = reduced_vector(reduced_index);
                ++reduced_index;
            }
        }
    }

    return expanded_matrix;
}

} } // namespace fem::mattools
