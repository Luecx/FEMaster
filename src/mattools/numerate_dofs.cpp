/**
 * @file numerate_dofs.cpp
 * @brief Implements contiguous numbering of active system degrees of freedom.
 *
 * The implementation is independent of a fixed nodal component count. It retains
 * the dimensions of the supplied activation mask and assigns compact algebraic
 * identifiers only to active entries, which allows structural node-by-six and
 * scalar thermal node-by-one systems to share the same utility.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 */

#include "numerate_dofs.h"

namespace fem { namespace mattools {

/**
 * Converts an arbitrary nodal activation mask into global system identifiers.
 *
 * The output preserves the exact row and component dimensions of `system_dofs`.
 * Active entries are numbered contiguously from zero in row-major traversal order
 * and inactive entries receive the sentinel `-1`. No structural six-component
 * assumption is made, so the same operation can enumerate scalar thermal systems.
 *
 * @param system_dofs Boolean node-by-component matrix identifying active unknowns.
 * @return System-index matrix with matching dimensions and contiguous active IDs.
 */
SystemDofIds numerate_dofs(const SystemDofs& system_dofs) {
    // Preserve the complete nodal/component layout of the activation mask so the
    // resulting mapping can be used directly by field reduction and assembly.
    SystemDofIds dof_ids(system_dofs.rows(), system_dofs.cols());
    int id_counter = 0;

    // Traverse in deterministic row-major order. Every active component receives
    // the next compact system index; inactive components retain the common -1
    // sentinel used throughout FEMaster's assembly utilities.
    for (Eigen::Index row = 0; row < system_dofs.rows(); ++row) {
        for (Eigen::Index col = 0; col < system_dofs.cols(); ++col) {
            dof_ids(row, col) = system_dofs(row, col) ? id_counter++ : -1;
        }
    }

    return dof_ids;
}

} } // namespace fem::mattools
