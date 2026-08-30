/**
 * @file numerate_dofs.cpp
 * @brief Implements contiguous numbering of active system DOFs.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 */

#include "numerate_dofs.h"

namespace fem { namespace mattools {

/**
 * @brief Generates SystemDofIds from an arbitrary node-by-component DOF mask.
 *
 * Active entries are numbered contiguously in row-major traversal order while
 * inactive entries receive -1. The returned matrix retains the exact dimensions
 * of the supplied mask, allowing structural node-by-six and scalar thermal
 * node-by-one systems to use the same numbering routine.
 *
 * @param system_dofs Boolean matrix indicating active nodal components.
 * @return Matrix of contiguous global DOF identifiers with matching dimensions.
 */
SystemDofIds numerate_dofs(const SystemDofs& system_dofs) {
    SystemDofIds dof_ids(system_dofs.rows(), system_dofs.cols());
    int id_counter = 0;

    for (Eigen::Index row = 0; row < system_dofs.rows(); ++row) {
        for (Eigen::Index col = 0; col < system_dofs.cols(); ++col) {
            dof_ids(row, col) = system_dofs(row, col) ? id_counter++ : -1;
        }
    }

    return dof_ids;
}

} } // namespace fem::mattools
