/**
 * @file reduce_mat_to_vec.h
 * @brief Provides functions to reduce a SystemDofIds matrix and NodeData into
 * a DynamicVector and to expand it back. Reduction ignores inactive DOFs, and
 * expansion fills the DynamicVector based on SystemDofIds.
 *
 * The `reduce_vec` function generates a reduced DynamicVector by extracting values
 * from the NodeData using the active DOFs indicated by the SystemDofIds matrix.
 *
 * The `expand_vec` function takes a DynamicVector and expands it to a NodeData
 * based on the active DOF indices in SystemDofIds.
 *
 * @author Created by Finn Eggers (c)
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"
#include <Eigen/Dense>

namespace fem { namespace mattools {

/**
 * @brief Reduces a NodeData into a DynamicVector by selecting active DOFs based
 * on the given SystemDofIds matrix.
 *
 * Inactive DOFs in SystemDofIds are marked by -1 and are ignored. The active
 * DOFs are collected into the output DynamicVector.
 *
 * @param dof_ids An Nx6 matrix that stores DOF indices, with -1 indicating inactive DOFs.
 * @param bc_matrix An Nx6 matrix that stores boundary condition values for each DOF.
 * @return DynamicVector The reduced vector containing only the active DOFs.
 */
DynamicVector reduce_mat_to_vec(const SystemDofIds& dof_ids, const NodeData& bc_matrix);

/**
 * @brief Expands a reduced DynamicVector back into a NodeData based on the
 * SystemDofIds matrix.
 *
 * Inactive DOFs (marked by -1) are skipped, while active DOFs are filled with
 * values from the DynamicVector.
 *
 * @param dof_ids An Nx6 matrix that stores DOF indices, with -1 indicating inactive DOFs.
 * @param reduced_vector The DynamicVector containing values for the active DOFs.
 * @return NodeData The expanded matrix containing the boundary conditions.
 */
NodeData expand_vec_to_mat(const SystemDofIds& dof_ids, const DynamicVector& reduced_vector);

} } // namespace fem::mattools
