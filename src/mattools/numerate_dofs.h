/******************************************************************************
* @file numerate_dofs.h
 * @brief numerate_dofs.h declares a function for numerating system DOF IDs
 * from a SystemDofs matrix in FEM simulations.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace mattools {

/******************************************************************************
 * @brief Generates SystemDofIds from a SystemDofs matrix by numerating
 * the degrees of freedom.
 *
 * This function takes a `SystemDofs` matrix (boolean values indicating DOFs)
 * and generates a `SystemDofIds` matrix with enumerated DOF IDs for each node.
 * DOF IDs start from 0 and are incremented for each DOF across the nodes.
 *
 * @param systemDofs The SystemDofs matrix to be converted.
 * @return SystemDofIds The matrix of numerated DOF IDs.
 ******************************************************************************/
SystemDofIds numerate_dofs(const SystemDofs& systemDofs);

} } // namespace fem::mattools
